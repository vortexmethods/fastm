/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.4    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/03/28     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: cuKernels.cuh                                                    |
| Info: Source code of BHgpu                                                  |
|                                                                             |
| This file is part of BHgpu.                                                 |
| BHcu is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BHgpu.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовки интерфейса с CUDA
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.4
\date 28 марта 2023 г.
*/

#ifndef CUKERNELS_CUH_
#define CUKERNELS_CUH_

#include "Params.h"

#ifdef __CUDACC__
#include <cuda.h>
#endif

#include <vector>

#include "Point2D.h"

#define THREADS1 32
#define THREADS2 512
#define THREADS3 64

#if __CUDA_ARCH__ >= 800
#define THREADS4 352
#else
#define THREADS4 192
#endif

#if __CUDA_ARCH__ >= 800
#define THREADS5 1024
#else
#define THREADS5 512
#endif

#define FACTOR1 6
#define FACTOR2 2
#define FACTOR3 4
#define FACTOR4 4
#define FACTOR5 2

namespace BHcu
{
    void CudaSelect(int dev);
    
    void setBlocks(int& blocks_);

    void* cudaNew(int n, size_t sizeType);

    void cudaDelete(void* cudaPtr);

    void cudaCopyVecToDevice(void* hostPtr, void* cudaPtr, size_t n, size_t typeSize);

    void cudaCopyVecFromDevice(void* cudaPtr, void* hostPtr, size_t n, size_t typeSize);


    void CudaTest(const char* msg);

    void KernelsOptimization();

    /******************************************************************************/
    /*** initialize memory ********************************************************/
    /******************************************************************************/

    float cuInitializationKernel(int* __restrict errd);


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/

    float cuBoundingBoxKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd, volatile int* __restrict childd,
        volatile int* __restrict massd,
        volatile real * __restrict momsd,
        volatile realPoint* __restrict posd,
        volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd);


	float McuBoundingBoxKernel(
		int nbodiesd,
		volatile realPoint* __restrict posd,
		volatile realPoint* __restrict Mposd,
		volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd);

	/******************************************************************************/
	/*** Morton codes *************************************************************/
	/******************************************************************************/

	float McuMortonCodesKernel(
        int nbodiesd,
        realPoint* __restrict posd,
        int* __restrict MmortonCodesKeyUnsortd, int* __restrict MmortonCodesIdxUnsortd,
        int* __restrict MmortonCodesKeyd, int* __restrict MmortonCodesIdxd,
        intPair* __restrict Mranged);


    /******************************************************************************/
    /*** Morton Internal nodes build **********************************************/
    /******************************************************************************/

    float McuMortonInternalNodesKernel(
        int nbodiesd,
        int* __restrict MmortonCodesKeyd,
        int* __restrict Mparentd,
        intPair* __restrict Mchildd,
        intPair* __restrict Mranged
    );

    /******************************************************************************/
    /*** Morton Internal nodes geometry calculation *******************************/
    /******************************************************************************/
    float McuMortonInternalCellsGeometryKernel(
        int nbodiesd,
        int* __restrict MmortonCodesKeyd,
        realPoint* __restrict Mposd,
        realPoint* __restrict Msized,
        intPair* __restrict Mranged,
        int* __restrict MlevelUnsortd, 
        int* __restrict MlevelSortd, 
        int* __restrict MindexUnsortd, 
        int* __restrict MindexSortd
    );

    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/
    float cuClearKernel23(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd,
        volatile int* __restrict massd,
        const real* __restrict gamd,
        volatile real* __restrict momsd);

    /******************************************************************************/
    /*** compute center of mass ***************************************************/
    /******************************************************************************/
    float cuSummarizationKernel2(
        const int nnodesd, const int nbodiesd,
        volatile int* __restrict countd, const int* __restrict childd,
        volatile int* __restrict massd,
        volatile real* __restrict momsd,
        volatile realPoint* __restrict posd,
        const int* __restrict cftl);
    /******************************************************************************/
    /*** sort bodies **************************************************************/
    /******************************************************************************/
    float cuSortKernel2(
        int nnodesd, int nbodiesd,
        volatile int* __restrict sortd, const int* __restrict countd,
        volatile int* __restrict startd, volatile int* __restrict childd);

    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/
    float cuForceCalculationKernel2(
        int nnodesd, int nbodiesd,
        int* __restrict errd,
        real itolsqd, real epssqd,
        const int* __restrict sortd, const int* __restrict childd,
        const real* __restrict momsd,
        const realPoint* __restrict posd,
        volatile realPoint* __restrict veld,
        volatile realPoint* __restrict Msized);

    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/
    float cuForceDirectCalculationKernel(
        int nnodesd, int nbodiesd,
        int* __restrict errd,
        real itolsqd, real epssqd,
        const int* __restrict sortd, const int* __restrict childd,
        const real * __restrict momsd,
        const realPoint* __restrict posd,
        volatile realPoint* __restrict veld);
}//namespace BHcu

#endif