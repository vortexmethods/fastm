/*--------------------------------*- BHcu -*-----------------*---------------*\
| #####   ##  ##                |                            | Version 1.0    |
| ##  ##  ##  ##   ####  ##  ## |  BHcu: Barnes-Hut method   | 2021/08/05     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: cuKernels.cuh                                                    |
| Info: Source code of BHcu                                                   |
|                                                                             |
| This file is part of BHcu.                                                  |
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
| along with BHcu.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовки интерфейса с CUDA
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
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
    struct moms
    {
        real gam;
#ifdef USE_DIP
        real foo;

#ifdef __CUDACC__
        real2 dip;
#else
        realPoint dip;
#endif

#ifdef USE_QUA

#ifdef __CUDACC__
        real2 qua;
#else
        realPoint qua;
#endif

       
#ifdef USE_OCT

#ifdef __CUDACC__
        real2 oct;
#else
        realPoint oct;
#endif

#ifdef USE_HEX
#ifdef __CUDACC__
        real2 hex;
#else
        realPoint hex;
#endif
#endif
#endif
#endif
#endif
    };

 
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
        volatile moms* __restrict momd,
        volatile realPoint* __restrict posd,
        volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd);


    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/

    float cuClearKernel1(int nnodesd, int nbodiesd, volatile int* __restrict childd);

    float cuTreeBuildingKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict errd, volatile int* __restrict childd,
        volatile realPoint* __restrict posd);

    float cuClearKernel23(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd,
        volatile int* __restrict massd,
        const real* __restrict gamd,
        volatile moms* __restrict momd);

    /******************************************************************************/
    /*** compute center of mass ***************************************************/
    /******************************************************************************/

    float cuSummarizationKernel(
        const int nnodesd, const int nbodiesd,
        volatile int* __restrict countd, const int* __restrict childd,
        volatile int* __restrict massd,
        volatile moms* __restrict momd,
        volatile realPoint* __restrict posd);


    /******************************************************************************/
    /*** sort bodies **************************************************************/
    /******************************************************************************/

    float cuSortKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict sortd, const int* __restrict countd,
        volatile int* __restrict startd, volatile int* __restrict childd);


    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/

    float cuForceCalculationKernel(
        int nnodesd, int nbodiesd,
        int* __restrict errd,
        real itolsqd, real epssqd,
        const int* __restrict sortd, const int* __restrict childd,
        const moms* __restrict momd,
        const realPoint* __restrict posd,
        volatile realPoint* __restrict veld);


    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/

    float cuForceDirectCalculationKernel(
        int nnodesd, int nbodiesd,
        int* __restrict errd,
        real itolsqd, real epssqd,
        const int* __restrict sortd, const int* __restrict childd,
        const moms* __restrict momd,
        const realPoint* __restrict posd,
        volatile realPoint* __restrict veld);
}//namespace BHcu

#endif