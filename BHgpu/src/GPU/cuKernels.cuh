/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.5    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/08/29     |
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
\version 1.5
\date 29 августа 2023 г.
*/

#ifndef CUKERNELS_CUH_
#define CUKERNELS_CUH_

#include "Params.h"

#ifdef __CUDACC__
#include <cuda.h>
#endif

#include <vector>

#include "Point2D.h"
#include "Vortex2D.h"

#define THREADS1 32
#define THREADS2 512
#define THREADS3 32
#if __CUDA_ARCH__ >= 800
#define THREADS4 352
#else
#define THREADS4 192
#endif

#if __CUDA_ARCH__ >= 800
#define THREADS5 1024
#else
#define THREADS5 1024
#endif

#define FACTOR1 6
#define FACTOR2 2
#define FACTOR3 4
#define FACTOR4 4
#define FACTOR5 1

namespace BHcu
{

    void setBinomCftConst(int* cft);
    
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

    float cuInitializationKernel();


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/

	float McuBoundingBoxKernel(
		int nbodiesd,
		const realVortex* __restrict vtxd,
		realPoint* __restrict Mposd,
		volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd);

	/******************************************************************************/
	/*** Morton codes *************************************************************/
	/******************************************************************************/

	float McuMortonCodesKernel(
        int nbodiesd,
        realVortex* __restrict vtxd,
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
        int nnodesd,
        int* __restrict MmortonCodesKeyd,
        realPoint* __restrict Mposd,
        realPoint* __restrict Msized,
        intPair* __restrict Mranged,
        int* __restrict MlevelUnsortd,
        int* __restrict MlevelSortd,
        int* __restrict MindexUnsortd,
        int* __restrict MindexSortd,
        int* __restrict MindexSortTd
    );
   


    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/
    float cuClearKernel2(
        int nnodesd, int nbodiesd,
        volatile int* __restrict massd,                
        volatile realPoint* __restrict momsd);

    /******************************************************************************/
    /*** compute center of mass ***************************************************/
    /******************************************************************************/
    float cuSummarizationKernel2(
        const int nnodesd, const int nbodiesd,
        const intPair* __restrict Mchildd,
        volatile int* __restrict massd,
        volatile realPoint* __restrict momsd,
        const realVortex* __restrict vtxd, const int* __restrict MmortonCodesIdxd,
        const realPoint* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd);

    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/
    float cuForceCalculationKernel2(
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const intPair* __restrict Mchildd,
        const realPoint* __restrict momsd,
        const realVortex* __restrict vtxd, const int* __restrict MmortonCodesIdxd,
        const realPoint* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
        volatile realPoint* __restrict veld,
        const realPoint* __restrict Msized);

    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/
    float cuForceDirectCalculationKernel(
        const int nnodesd, const int nbodiesd,
        const real epssqd,
        const realVortex* __restrict vtxd,        
        volatile realPoint* __restrict veld);
}//namespace BHcu

#endif