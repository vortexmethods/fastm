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
| File name: cuKernels.cu                                                     |
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
\brief Реализация CUDA-ядер
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#include "cuKernels.cuh"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <cuda.h>
#include "operations.cuh"
#include "cuSort.cuh"

#include "Point2D.h"

#define WARPSIZE 32
#define MAXDEPTH 28
#define BLOCKD 32

#define codeLength 14
#define twoPowCodeLength (1 << codeLength)

namespace BHcu
{
    int blocks;
__device__ volatile int bottomd;
__device__ unsigned int blkcntd;
__device__ volatile real radiusd;

void setBlocks(int& blocks_)
{
     blocks_ = blocks;
}

void CudaSelect(int dev)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        fprintf(stderr, "There is no device supporting CUDA\n");
        exit(-1);
    }

    if ((dev < 0) || (deviceCount <= dev)) {
        fprintf(stderr, "There is no device %d\n", dev);
        exit(-1);
    }
    cudaSetDevice(dev);


    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if ((deviceProp.major == 9999) && (deviceProp.minor == 9999)) {
        fprintf(stderr, "There is no CUDA capable device\n");
        exit(-1);
    }
    if (deviceProp.major < 2) {
        fprintf(stderr, "Need at least compute capability 2.0\n");
        exit(-1);
    }
    if (deviceProp.warpSize != WARPSIZE) {
        fprintf(stderr, "Warp size must be %d\n", deviceProp.warpSize);
        exit(-1);
    }

    blocks = deviceProp.multiProcessorCount;
    
    if ((WARPSIZE <= 0) || (WARPSIZE & (WARPSIZE - 1) != 0)) {
        fprintf(stderr, "Warp size must be greater than zero and a power of two\n");
        exit(-1);
    }
    if (MAXDEPTH > WARPSIZE) {
        fprintf(stderr, "MAXDEPTH must be less than or equal to WARPSIZE\n");
        exit(-1);
    }
    if ((THREADS1 <= 0) || (THREADS1 & (THREADS1 - 1) != 0)) {
        fprintf(stderr, "THREADS1 must be greater than zero and a power of two\n");
        exit(-1);
    }

    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, 0);

    int fact = 1024;
    int driverVersion, runtimeVersion;

    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);

    printf("\n");
    printf("                          GPU Device Properties                         \n");
    printf("------------------------------------------------------------------------\n");
    printf("Name:                                  %s\n", properties.name); 
    printf("CUDA driver/runtime version:           %d.%d/%d.%d\n", driverVersion / 1000, (driverVersion % 100) / 10, runtimeVersion / 1000, (runtimeVersion % 100) / 10);
    printf("CUDA compute capabilitiy:              %d.%d\n", properties.major, properties.minor);
    printf("Number of multiprocessors:             %d\n", properties.multiProcessorCount);
    if (printFullCUDAinfo)
    {
        printf("GPU clock rate:                        %d (MHz)\n", properties.clockRate / fact);
        printf("Memory clock rate:                     %d (MHz)\n", properties.memoryClockRate / fact);
        printf("Memory bus width:                      %d-bit\n", properties.memoryBusWidth);
        printf("Theoretical memory bandwidth:          %d (GB/s)\n", (properties.memoryClockRate / fact * (properties.memoryBusWidth / 8) * 2) / fact);
        printf("Device global memory:                  %d (MB)\n", (int)(properties.totalGlobalMem / (fact * fact)));
        printf("Shared memory per block:               %d (KB)\n", (int)(properties.sharedMemPerBlock / fact));
        printf("Constant memory:                       %d (KB)\n", (int)(properties.totalConstMem / fact));
        printf("Maximum number of threads per block:   %d\n", properties.maxThreadsPerBlock);
        printf("Maximum thread dimension:              [%d, %d, %d]\n", properties.maxThreadsDim[0], properties.maxThreadsDim[1], properties.maxThreadsDim[2]);
        printf("Maximum grid size:                     [%d, %d, %d]\n", properties.maxGridSize[0], properties.maxGridSize[1], properties.maxGridSize[2]);
    }
    printf("------------------------------------------------------------------------\n");  
}


void cudaDelete(void* cudaPtr)
{
	cudaFree(cudaPtr);
}


void* cudaNew(int n, size_t sizeType)
{
	void* cudaPtr;
	cudaMalloc(&cudaPtr, sizeType * n);
	CudaTest("couldn't allocate device memory");

	return cudaPtr;
}

void cudaCopyVecToDevice(void* hostPtr, void* cudaPtr, size_t n, size_t typeSize)
{
	cudaMemcpy(cudaPtr, hostPtr, typeSize * n, cudaMemcpyHostToDevice);
	CudaTest("couldn't copy data from host to device");
}

void cudaCopyVecFromDevice(void* cudaPtr, void* hostPtr, size_t n, size_t typeSize)
{
	cudaMemcpy(hostPtr, cudaPtr, typeSize * n, cudaMemcpyDeviceToHost);
	CudaTest("couldn't copy data from device to host");
}


//////////////////
/// Error TEST
//////////////////


void CudaTest(const char* msg)
{
    cudaError_t e;

    //cudaThreadSynchronize();
    cudaDeviceSynchronize();
    if (cudaSuccess != (e = cudaGetLastError())) {
        fprintf(stderr, "%s: %d\n", msg, e);
        fprintf(stderr, "%s\n", cudaGetErrorString(e));
        exit(-1);
    }
}


//////////////////
/// CUDA Kernels
//////////////////

__constant__ int binomCft[order * (order + 1)];

void setBinomCftConst(int* cft)
{
    cudaMemcpyToSymbol(binomCft, cft, order * (order + 1) * sizeof(int));
}



/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void InitializationKernel()
{
    blkcntd = 0;
}


/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS1, FACTOR1)
void MBoundingBoxKernel(
    const int nbodiesd, 
    const real3* __restrict vtxd, 
    real2* __restrict Mposd, 
    volatile real2* __restrict maxrd, 
    volatile real2* __restrict minrd)
{
    register int i, j, k, inc;
    register real2 val;
    register real2 minr, maxr;
    __shared__ volatile real2 sminr[THREADS1], smaxr[THREADS1];

    // initialize with valid data (in case #bodies < #threads)
    minr.x = maxr.x = vtxd[0].x;
    minr.y = maxr.y = vtxd[0].y;

    // scan all bodies
    i = threadIdx.x;
    inc = THREADS1 * gridDim.x;
    for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) {
        val.x = vtxd[j].x;
        val.y = vtxd[j].y;

        minr.x = realmin(minr.x, val.x);
        maxr.x = realmax(maxr.x, val.x);

        minr.y = realmin(minr.y, val.y);
        maxr.y = realmax(maxr.y, val.y);
    }

    // reduction in shared memory
    sminr[i].x = minr.x;
    smaxr[i].x = maxr.x;
    sminr[i].y = minr.y;
    smaxr[i].y = maxr.y;

    for (j = THREADS1 / 2; j > 0; j /= 2) {
        __syncthreads();
        if (i < j) {
            k = i + j;
            sminr[i].x = minr.x = realmin(minr.x, sminr[k].x);
            smaxr[i].x = maxr.x = realmax(maxr.x, smaxr[k].x);
            sminr[i].y = minr.y = realmin(minr.y, sminr[k].y);
            smaxr[i].y = maxr.y = realmax(maxr.y, smaxr[k].y);
        }
    }

    // write block result to global memory
    if (i == 0) {
        k = blockIdx.x;
        minrd[k].x = minr.x;
        maxrd[k].x = maxr.x;
        minrd[k].y = minr.y;
        maxrd[k].y = maxr.y;

        

        __threadfence();

        inc = gridDim.x - 1;
        if (inc == atomicInc(&blkcntd, inc)) {
            // I'm the last block, so combine all block results
            for (j = 0; j <= inc; j++) {
                minr.x = realmin(minr.x, minrd[j].x);
                maxr.x = realmax(maxr.x, maxrd[j].x);
                minr.y = realmin(minr.y, minrd[j].y);
                maxr.y = realmax(maxr.y, maxrd[j].y);
            }

            // compute 'radius'
            radiusd = realmax(maxr.x - minr.x, maxr.y - minr.y) / 2;

            // create root node
            Mposd[0].x = (minr.x + maxr.x) / 2;
            Mposd[0].y = (minr.y + maxr.y) / 2;            
        }
    }
}

/******************************************************************************/
/*** Morton codes *************************************************************/
/******************************************************************************/

__global__
void MMortonCodesKernel (
    const int nbodies, 
    const real3* __restrict vtxd, 
    int* __restrict MmortonCodesKeyUnsortd, 
    int* __restrict MmortonCodesIdxUnsortd)
{
	int bdy = blockDim.x * blockIdx.x + threadIdx.x;

	if (bdy < nbodies)
	{
		real x = twoPowCodeLength * vtxd[bdy].x;
		real y = twoPowCodeLength * vtxd[bdy].y;

		unsigned int xx = MExpandBits((unsigned int)x);
		unsigned int yy = MExpandBits((unsigned int)y);
		MmortonCodesKeyUnsortd[bdy] = yy | (xx << 1);
		MmortonCodesIdxUnsortd[bdy] = bdy;
	}
}


/******************************************************************************/
/*** Morton Internal nodes tree build *****************************************/
/******************************************************************************/
__global__
void MMortonInternalNodesKernel(
    const int nbodies, 
    const int* __restrict MmortonCodesKeyd, 
    int* __restrict Mparentd, 
    int2* __restrict Mchildd, 
    int2* __restrict Mranged)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (i < nbodies - 1)
    {
        int d = sign(Delta(i, i + 1, nbodies, MmortonCodesKeyd) - Delta(i, i - 1, nbodies, MmortonCodesKeyd));
        int delta_min = Delta(i, i - d, nbodies, MmortonCodesKeyd);

        int Lmax = 2;
        while (Delta(i, i + Lmax * d, nbodies, MmortonCodesKeyd) > delta_min)
            Lmax *= 2;

        int L = 0;
        for (int t = (Lmax >> 1); t >= 1; t >>= 1)
            if (Delta(i, i + (L + t) * d, nbodies, MmortonCodesKeyd) > delta_min)
                L += t;

        int j = i + L * d;

        int delta_node = Delta(i, j, nbodies, MmortonCodesKeyd);

        int s = 0;
        for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
        {
            int dl = Delta(i, i + (s + t) * d, nbodies, MmortonCodesKeyd);
            if (dl > delta_node)
                s += t;
        }//for p


        int gamma = i + s * d +   d * (d < 0);   //последнее слагаемое = std::min(d, 0);

        int Mmin = min(i, j);
        int Mmax = max(i, j);
        
        const int& left = gamma;
        const int& right = gamma + 1;

        // Левый потомок - лист или внутренний узел
        int childLeft = Mchildd[i].x = (Mmin == gamma) * nbodies + left;
        
        Mranged[childLeft].x = Mmin;
        Mranged[childLeft].y = gamma;
        Mparentd[childLeft] = i;

        // Правый потомок - лист или внутренний узел
        int childRight = Mchildd[i].y = (Mmax == gamma + 1) * nbodies + right;

        Mranged[childRight].x = gamma+1;
        Mranged[childRight].y = Mmax;
        Mparentd[childRight] = i;
    }
}

/******************************************************************************/
/*** Morton Internal nodes geometry calculation *******************************/
/******************************************************************************/
__global__
void MMortonInternalCellsGeometryKernel(
    const int nbodies,
    const int* __restrict MmortonCodesKeyd,
    real2* __restrict Mposd,
    real2* __restrict Msized,
    const int2* __restrict Mranged,
    int* __restrict MlevelUnsortd,
    int* __restrict MindexUnsortd
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;

    if (cell < nbodies - 1)
    {
        int prLength = min(Delta(Mranged[cell].x, Mranged[cell].y, nbodies, MmortonCodesKeyd), 2 * codeLength);
        unsigned int pr = (MmortonCodesKeyd[Mranged[cell].x] >> (2 * codeLength - prLength));
               
        prLength -= min(Delta(Mranged[0].x, Mranged[0].y, nbodies, MmortonCodesKeyd), 2 * codeLength);
        
        real2 sz;
        sz.x = 1 / (real)(1 << ceilhalf(prLength));
        sz.y = 1 / (real)(1 << (prLength / 2));

        real2 pos;
        pos.x = sz.x / 2;
        pos.y = sz.y / 2;

        int xint = MShrinkBits(pr);       
        int yint = MShrinkBits(pr >> 1);

        real addX = xint * sz.x;
        real addY = yint * sz.y;

        if (prLength & 1)
        {
            pos.x += addX;
            pos.y += addY;
        }
        else
        {
            pos.y += addX;
            pos.x += addY;
        }			          

        Mposd[cell] = pos;   
        Msized[cell] = sz;    

                
        MlevelUnsortd[cell] = prLength;
        MindexUnsortd[cell] = cell;
    }

}//MMortonInternalCellsGeometryKernel(...)


/******************************************************************************/
/*** permutation list transposition *******************************************/
/******************************************************************************/
__global__
void MTransposeIndexKernel(
    const int nbodiesd, const int nnodesd,
    const int* __restrict MindexSortd, 
    int* __restrict MindexSortTd)
{
    register const int cell = blockDim.x * blockIdx.x + threadIdx.x;
    register const int newcell = MindexSortd[cell];

    if (cell < nbodiesd - 1)
        MindexSortTd[newcell] = cell;       

}//MTransposeIndexKernel



/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(1024, 1)
void ClearKernel2(
    const int nnodesd, const int nbodiesd, 
    volatile int* __restrict massd)
{
    register int k, inc, bottom;

    bottom = nnodesd - (nbodiesd - 1); //bottomd;
    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    if (k < bottom) k += inc;

// 0 1 ... (nb-1)  (nb+0) ... (nb+nb-2)
// --------------  --------------------
//     bodies              cells

    // iterate over all cells assigned to thread
    while (k < nnodesd) 
	{
        massd[nnodesd - 1 - k] = -1;
        k += inc;
    }
}



/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/
#include "ShiftKernels/IncludeKer.cu"

/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS5, FACTOR5)
void ForceCalculationKernel2(
    const int nnodesd, const int nbodiesd,
    const real itolsqd, const real epssqd,
    const int2* __restrict Mchildd,
    const real2* __restrict momsd,
    const real3* __restrict vtxd,    
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    real2* __restrict veld,  //veld - без volatile
    const real2* __restrict Msized)

{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, v, dr, ps;
    register real r2;
    register const real2* mom;
    //register real2 mom[order];  
    //register real2 mom0, mom1, mom2, mom3, mom4, mom5, mom6, mom7, mom8, mom9, mom10, mom11, mom12, mom13;

    register real2 th;
    

    __shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE], node[MAXDEPTH * THREADS5 / WARPSIZE];



    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;
    //diff = threadIdx.x - sbase;

    __syncthreads();
    __threadfence_block();

    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
    {
        const int indexInParticles = MmortonCodesIdxd[k];
        p = real2{ vtxd[indexInParticles].x, vtxd[indexInParticles].y };

        v.x = 0;
        v.y = 0;

        // initialize iteration stack, i.e., push root node onto stack
        depth = j;
        if (sbase == threadIdx.x)
        {
            pos[j] = 0;
            node[j] = nnodesd - 1;
        }

        do
        {
            // stack is not empty
            pd = pos[depth];
            nd = node[depth];

            register int2 chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

			register real gm;
			register real2 sumSide2;
			bool isVortex;

            while (pd < 2)
            {
                // node on top of stack has more children to process

                // load child pointer
                //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
				int chd = pd * chBoth.y + (1-pd) * chBoth.x;
				++pd;
				
				isVortex = (chd >= nbodiesd);
				
				if (isVortex)
				{
					n = chd - nbodiesd;
					ps = real2{ vtxd[MmortonCodesIdxd[n]].x, vtxd[MmortonCodesIdxd[n]].y };
					gm = vtxd[MmortonCodesIdxd[n]].z;
					sumSide2 = real2{ (real)0, (real)0 };
				}
				else
				{
					register const int srtT = MindexSortTd[chd];
					n = (nnodesd - 1) - srtT;
					ps = Mposd[chd];
					mom = momsd + (srtT * order);
					gm = mom[0].x;
					sumSide2 = Msized[chd];
				}

				dr = p - ps;

				//for (i = 0; i < order; ++i)
				//    mom[i] = momsd[n * order + i];

				//mom0 = momsd[n * order];
				//mom1 = momsd[n * order + 1];
				//mom2 = momsd[n * order + 2];
				//mom3 = momsd[n * order + 3];
				//mom4 = momsd[n * order + 4];
				//mom5 = momsd[n * order + 5];
				//mom6 = momsd[n * order + 6];
				//mom7 = momsd[n * order + 7];
				//mom8 = momsd[n * order + 8];
				//mom9 = momsd[n * order + 9];
				//mom10 = momsd[n * order + 10];
				//mom11 = momsd[n * order + 11];
				//mom12 = momsd[n * order + 12];
				//mom13 = momsd[n * order + 13];

				r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared

				// check if all threads agree that cell is far enough away (or is a body)
				if (isVortex || __all_sync(0xffffffff, ((sumSide2.x+sumSide2.y)*(sumSide2.x+sumSide2.y) + epssqd) * itolsqd < r2))
				{  
#ifdef CALCinDOUBLE
					real f = gm / realmax(r2, epssqd);
#else
					real f = fdividef(gm, realmax(r2,epssqd));
#endif
					v += f * dr;

					if ((!isVortex) && (order > 1))
					{
#ifdef CALCinDOUBLE
						real2 cftr = (r2 ? (1.0 / r2) : (real)0) * dr;
#else
						real2 cftr = (r2 ? fdividef(1.0f, r2) : 0.0f) * dr;
#endif
						th = cftr;
						
						for (int s = 1; s < order; ++s)
						{
							th = multz(th, cftr);
#ifdef CALCinFLOAT                                    
							if (isinf(th.x) || isinf(th.y))
							{
								//printf("s = %d\n", s);
								break;
							}
#endif
							v += multzA(th, mom[s]);
						}
						
						//th = multz(th, cftr);
						//v += multzA(th, mom1);

						//th = multz(th, cftr);
						//v += multzA(th, mom2);

						//th = multz(th, cftr);
						//v += multzA(th, mom3);

						//th = multz(th, cftr);
						//v += multzA(th, mom4);

						//th = multz(th, cftr);
						//v += multzA(th, mom5);

						//th = multz(th, cftr);
						//v += multzA(th, mom6);

						//th = multz(th, cftr);
						//v += multzA(th, mom7);

						//th = multz(th, cftr);
						//v += multzA(th, mom8);

						//th = multz(th, cftr);
						//v += multzA(th, mom9);

						//th = multz(th, cftr);
						//v += multzA(th, mom10);

						//th = multz(th, cftr);
						//v += multzA(th, mom11);

						//th = multz(th, cftr);
						//v += multzA(th, mom12);

						//th = multz(th, cftr);
						//v += multzA(th, mom13);
					}
				}
				else
				{
					// push cell onto stack
					if (sbase == threadIdx.x)
					{  // maybe don't push and inc if last child
						pos[depth] = pd;
						node[depth] = nd;
					}
					depth++;
					pd = 0;
					nd = n;

					chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
				}
                
            }
            depth--;  // done with this level
        } while (depth >= j);

        // update velocity

        real2 result = real2{ -v.y, v.x };
        veld[indexInParticles] = result;        
    }
} 


/******************************************************************************/
/*** compute force (direct) ***************************************************/
/******************************************************************************/


__global__
//__launch_bounds__(THREADSD, FACTORD)
void ForceDirectCalculationKernel(
    const int nnodesd, const int nbodiesd,
    const real epssqd,
    const real3* __restrict vtxd,    
    real2* __restrict veld)
{
    __shared__ real3 shvtx[BLOCKD];    

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    real2 pt;
    pt.x = vtxd[i].x;
    pt.y = vtxd[i].y;

    real2 vel;
    vel.x = vel.y = 0;

    real2 dr;
    real dr2, izn;

    //vortices
    for (size_t j = 0; j < nbodiesd; j += BLOCKD)
    {
        shvtx[threadIdx.x] = vtxd[j + threadIdx.x];               

        __syncthreads();

        for (size_t q = 0; q < BLOCKD; ++q)
        {
            if (j + q < nbodiesd)
            {
                dr.x = pt.x - shvtx[q].x;
                dr.y = pt.y - shvtx[q].y;
                dr2 = dr.x * dr.x + dr.y * dr.y;

                izn = shvtx[q].z / realmax(dr2, epssqd);// / CUboundDenom(dr2, eps2); //РЎРіР»Р°Р¶РёРІР°С‚СЊ РЅР°РґРѕ!!!

                vel.x -= dr.y * izn;
                vel.y += dr.x * izn;

            }
        }
        __syncthreads();
    }

    if (i < nbodiesd)
    {
        veld[i].x = vel.x;// * iDPI;
        veld[i].y = vel.y;// * iDPI;
    }
    //*/
}


/******************************************************************************/


void KernelsOptimization()
{
    // set L1/shared memory configuration

    cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value
        
    cudaFuncSetCacheConfig(ForceDirectCalculationKernel, cudaFuncCachePreferEqual); //d
    cudaGetLastError();  // reset error value
        
    cudaFuncSetCacheConfig(SummarizationKernel2_14, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_16, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ForceCalculationKernel2, cudaFuncCachePreferL1); //d
    cudaGetLastError();  // reset error value

}


/******************************************************************************/





//////////////////
/// Wrappers
//////////////////



    /******************************************************************************/
    /*** initialize memory ********************************************************/
    /******************************************************************************/

    float cuInitializationKernel()
    {
        //fprintf(stderr, "IKKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);

        cudaEventRecord(start, 0);
        InitializationKernel<<<1, 1>>> ();
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        CudaTest("kernel 0 launch failed");
        
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/
	float McuBoundingBoxKernel(
		int nbodiesd,
		const realVortex* __restrict vtxd,
		realPoint* __restrict Mposd,
		volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd)
	{
		cudaEvent_t start, stop;
		float time;

		cudaEventCreate(&start);  cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		MBoundingBoxKernel<<<blocks * FACTOR1, THREADS1>>> (nbodiesd, (real3*)vtxd, (real2*)Mposd, (real2*)maxrd, (real2*)minrd);
		cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

		CudaTest("Mkernel 1 launch failed");

		cudaEventDestroy(start);  cudaEventDestroy(stop);
		return time;
	}

	/******************************************************************************/
	/*** Morton codes *************************************************************/
	/******************************************************************************/

    float McuMortonCodesKernel(
        int nbodiesd,
        realVortex* __restrict vtxd,
        int* __restrict MmortonCodesKeyUnsortd, int* __restrict MmortonCodesIdxUnsortd,
        int* __restrict MmortonCodesKeyd, int* __restrict MmortonCodesIdxd,
        intPair* __restrict Mranged)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = (nbodiesd + 31) / 32;
        dim3 Mthreads = 32;

        MMortonCodesKernel << <Mblocks, Mthreads >> > (nbodiesd, (real3*)vtxd, MmortonCodesKeyUnsortd, MmortonCodesIdxUnsortd);


        ///RadixSort

        RadixSortFromCUB(
            MmortonCodesKeyUnsortd, MmortonCodesKeyd, \
            MmortonCodesIdxUnsortd, MmortonCodesIdxd, \
            nbodiesd, 0, 2 * codeLength);


        //Заполнение нулевой ячейки (диапазон для корня дерева)
        int totalRange[2] = { 0, nbodiesd - 1 };
        cudaCopyVecToDevice(totalRange, Mranged, 2, sizeof(int));

		cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

		CudaTest("Mkernel 1 launch failed");

		cudaEventDestroy(start);  cudaEventDestroy(stop);

        

		return time;
	}

    /******************************************************************************/
    /*** Morton Internal nodes build **********************************************/
    /******************************************************************************/

    float McuMortonInternalNodesKernel(
        int nbodiesd,
        int* __restrict MmortonCodesKeyd, 
        int* __restrict Mparentd,
        intPair* __restrict Mchildd,
        intPair* __restrict Mranged
    )
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((nbodiesd - 1) + 31) / 32;
        dim3 Mthreads = 32;

        MMortonInternalNodesKernel<<<Mblocks, Mthreads>>> (nbodiesd, MmortonCodesKeyd, Mparentd, (int2*)Mchildd, (int2*)Mranged);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;

    }



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
    )
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((nbodiesd) + 31) / 32;
        dim3 Mthreads = 32;

		MMortonInternalCellsGeometryKernel<<<Mblocks, Mthreads>>>(nbodiesd, MmortonCodesKeyd, (real2*)Mposd, (real2*)Msized, (int2*)Mranged, MlevelUnsortd, MindexUnsortd);

        RadixSortFromCUB( \
            MlevelUnsortd, MlevelSortd, \
            MindexUnsortd, MindexSortd, \
            nbodiesd-1, 0, 2 * codeLength);

        MTransposeIndexKernel << <Mblocks, Mthreads >> > (nbodiesd, nnodesd, MindexSortd, MindexSortTd);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }



	 	 

    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/   

    float cuClearKernel2(
        int nnodesd, int nbodiesd,
        volatile int* __restrict massd,        
        volatile realPoint* __restrict momsd)
    {
        //fprintf(stderr, "CxKernel\n");
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        cudaMemset((void*)momsd, 0, (nbodiesd-1) * order * sizeof(realPoint));

        ClearKernel2 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, massd);
        
	       	

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel clear2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    /******************************************************************************/
    /*** compute multipole moments for all the cells ******************************/
    /******************************************************************************/
    float cuSummarizationKernel2(
        const int nnodesd, const int nbodiesd,
        const intPair* __restrict Mchildd,
        volatile int* __restrict massd,
        volatile realPoint* __restrict momsd,
        const realVortex* __restrict vtxd, const int* __restrict MmortonCodesIdxd,
        const realPoint* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd)
    {
        //fprintf(stderr, "SKKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        
#include "ShiftKernels/SwitchKer.cu"

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }

   

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
        const realPoint* __restrict Msized)
    {
        //fprintf(stderr, "FCKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ForceCalculationKernel2 << <blocks * FACTOR5, THREADS5 >> > (
            nnodesd, nbodiesd, itolsqd, epssqd, (int2*)Mchildd, (real2*)momsd,
            (real3*)vtxd, MmortonCodesIdxd,
            (real2*)Mposd, MindexSortd, MindexSortTd,
            (real2*)veld, (real2*)Msized);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }

    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/

    float cuForceDirectCalculationKernel(
        int nnodesd, int nbodiesd,
        real epssqd,
        const realVortex* __restrict vtxd,
        volatile realPoint* __restrict veld)
    {
        //fprintf(stderr, "DFKernel\n");
        
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);
        
        ForceDirectCalculationKernel<<<(nbodiesd + BLOCKD - 1) / BLOCKD, BLOCKD>>> (nnodesd, nbodiesd, epssqd, (real3*)vtxd, (real2*)veld);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);
        
        CudaTest("kernel direct launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);
        
        return timeD;
    }

}//namespace BHcu