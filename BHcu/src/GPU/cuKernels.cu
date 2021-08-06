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
| File name: cuKernels.cu                                                     |
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
\brief Реализация CUDA-ядер
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
*/



#include "cuKernels.cuh"

#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>


#include "operations.cuh"
#include "Point2D.h"

#define WARPSIZE 32
#define MAXDEPTH 32
#define BLOCKD 32


namespace BHcu
{

int blocks;
__device__ volatile int stepd, bottomd, maxdepthd;
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


/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void InitializationKernel(int* __restrict errd)
{
    *errd = 0;
    stepd = -1;
    maxdepthd = 1;
    blkcntd = 0;
}


/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS1, FACTOR1)
void BoundingBoxKernel(int nnodesd, int nbodiesd, volatile int* __restrict startd, volatile int* __restrict childd, volatile int* __restrict massd, volatile moms* __restrict momd, volatile real2* __restrict posd, volatile real2* __restrict maxrd, volatile real2* __restrict minrd)
{
    register int i, j, k, inc;
    register real2 val;
    register real2 minr, maxr;
    __shared__ volatile real2 sminr[THREADS1], smaxr[THREADS1];

    // initialize with valid data (in case #bodies < #threads)
    minr.x = maxr.x = posd[0].x;
    minr.y = maxr.y = posd[0].y;

    // scan all bodies
    i = threadIdx.x;
    inc = THREADS1 * gridDim.x;
    for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) {
        val.x = posd[j].x;
        val.y = posd[j].y;

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
            k = nnodesd;
            bottomd = k;

            massd[k] = -1;
            momd[k].gam = 0;
#ifdef USE_DIP
            momd[k].dip.x = 0;
            momd[k].dip.y = 0;
#ifdef USE_QUA
            momd[k].qua.x = 0;
            momd[k].qua.y = 0;
#ifdef USE_OCT
            momd[k].oct.x = 0;
            momd[k].oct.y = 0;
#ifdef USE_HEX            
            momd[k].hex.x = 0;
            momd[k].hex.y = 0;
#endif
#endif
#endif
#endif
            startd[k] = 0;
            posd[k].x = (minr.x + maxr.x) / 2;
            posd[k].y = (minr.y + maxr.y) / 2;
            k *= 4;
            for (i = 0; i < 4; i++) childd[k + i] = -1;

            stepd++;
        }
    }
}





/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(1024, 1)
void ClearKernel1(int nnodesd, int nbodiesd, volatile int* __restrict childd)
{
    register int k, inc, top, bottom;

    top = 4 * nnodesd;
    bottom = 4 * nbodiesd;
    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    if (k < bottom) k += inc;

    // iterate over all cells assigned to thread
    while (k < top) {
        childd[k] = -1;
        k += inc;
    }
}


__global__
__launch_bounds__(THREADS2, FACTOR2)
void TreeBuildingKernel(int nnodesd, int nbodiesd, volatile int* __restrict errd, volatile int* __restrict childd, volatile real2* __restrict posd)
{
    register int i, j, depth, localmaxdepth, skip, inc;
    register real x, y, r;
    register real2 p;
    register real2 dr;
    register int ch, n, cell, locked, patch;
    register real radius;
    register real2 root;

    // cache root data
    radius = radiusd;

	//if ((threadIdx.x == 0) && (blockIdx.x == 0))
	//	printf("radius = %f\n", radius);
	//	printf("pod[10].g = %f, %f\n", posd[10].x, posd[10].y);
	//	printf("nn = %d, %d\n", nnodesd, nbodiesd);

    root.x = posd[nnodesd].x;
    root.y = posd[nnodesd].y;

    localmaxdepth = 1;
    skip = 1;
    inc = blockDim.x * gridDim.x;
    i = threadIdx.x + blockIdx.x * blockDim.x;

    // iterate over all bodies assigned to thread
    while (i < nbodiesd) {
        p.x = posd[i].x;
        p.y = posd[i].y;

        if (skip != 0) {
            // new body, so start traversing at root
            skip = 0;
            p.x = posd[i].x;
            p.y = posd[i].y;

            n = nnodesd;
            depth = 1;
            r = radius / 2;
            dr.x = dr.y = -r;
            j = 0;
            // determine which child to follow
            if (root.x < p.x) { j = 1; dr.x = r; }
            if (root.y < p.y) { j |= 2; dr.y = r; }
            x = root.x + dr.x;
            y = root.y + dr.y;
        }

        // follow path to leaf cell
        ch = childd[n * 4 + j];
        while (ch >= nbodiesd) {
            n = ch;
            depth++;
            r /= 2;
            dr.x = dr.y = -r;
            j = 0;
            // determine which child to follow
            if (x < p.x) { j = 1; dr.x = r; }
            if (y < p.y) { j |= 2; dr.y = r; }
            x += dr.x;
            y += dr.y;
            ch = childd[n * 4 + j];
        }

        if (ch != -2) {  // skip if child pointer is locked and try again later
            locked = n * 4 + j;
            if (ch == -1) {
                if (-1 == atomicCAS((int*)&childd[locked], -1, i)) {  // if null, just insert the new body
                    localmaxdepth = max(depth, localmaxdepth);
                    i += inc;  // move on to next body
                    skip = 1;
                }
            }
            else {  // there already is a body in this position
                if (ch == atomicCAS((int*)&childd[locked], ch, -2)) {  // try to lock
                    patch = -1;
                    // create new cell(s) and insert the old and new body
                    do {
                        depth++;

                        cell = atomicSub((int*)&bottomd, 1) - 1;
                        if (cell <= nbodiesd) {
                            *errd = 1;
							//printf("!!!");
                            bottomd = nnodesd;
                        }

                        if (patch != -1) {
                            childd[n * 4 + j] = cell;
                        }
                        patch = max(patch, cell);

                        j = 0;
                        if (x < posd[ch].x) j = 1;
                        if (y < posd[ch].y) j |= 2;
                        childd[cell * 4 + j] = ch;

                        posd[cell].x = x;
                        posd[cell].y = y;


                        n = cell;
                        r /= 2;
                        dr.x = dr.y = -r;
                        j = 0;
                        if (x < p.x) { j = 1; dr.x = r; }
                        if (y < p.y) { j |= 2; dr.y = r; }
                        x += dr.x;
                        y += dr.y;

                        ch = childd[n * 4 + j];
                        // repeat until the two bodies are different children
                    } while (ch >= 0);
                    childd[n * 4 + j] = i;

                    localmaxdepth = max(depth, localmaxdepth);
                    i += inc;  // move on to next body
                    skip = 2;
                }
            }
        }
        __syncthreads();
        __threadfence();

        if (skip == 2) {
            childd[locked] = patch;
        }
    }
    // record maximum tree depth
    atomicMax((int*)&maxdepthd, localmaxdepth);
}


__global__
__launch_bounds__(1024, 1)
void ClearKernel2(int nnodesd, volatile int* __restrict startd, volatile int* __restrict massd)
{
    register int k, inc, bottom;

    bottom = bottomd;
    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    if (k < bottom) k += inc;

    // iterate over all cells assigned to thread
    while (k < nnodesd) {
        massd[k] = -1;
        startd[k] = -1;
        k += inc;
    }
}


__global__
__launch_bounds__(1024, 1)
void ClearKernel3(int nnodesd, int nbodies, const real* __restrict gamd, volatile moms* __restrict momd)
{
    register int k, inc;

    inc = blockDim.x * gridDim.x;
    k = threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    
    // iterate over all cells assigned to thread
    while (k < nnodesd) {        
        momd[k].gam = (k < nbodies) ? gamd[k] : 0;
#ifdef USE_DIP
        momd[k].dip.x = 0;
        momd[k].dip.y = 0;
#ifdef USE_QUA
        momd[k].qua.x = 0;
        momd[k].qua.y = 0;
#ifdef USE_OCT
        momd[k].oct.x = 0;
        momd[k].oct.y = 0;
#ifdef USE_HEX
        momd[k].hex.x = 0;
        momd[k].hex.y = 0;
#endif
#endif
#endif
#endif        
        k += inc;
    }
}





/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS3, FACTOR3)
void SummarizationKernel(const int nnodesd, const int nbodiesd, volatile int* __restrict countd, const int* __restrict childd, volatile int* __restrict massd, volatile moms* __restrict momd, volatile real2* __restrict posd)
{
    register int i, j, k, ch, inc, cnt, bottom, flag;
    register real g, cg;
#ifdef USE_DIP
    register real2 cen, dr;
    register real2 d, cd, mh;
#ifdef USE_QUA
    register real2 q, cq, dh;
#ifdef USE_OCT
    register real2 o, co, qh;
#ifdef USE_HEX
    register real2 h, chex, oh;
#endif
#endif
#endif
#endif
    register int m, cm;
    __shared__ int child[THREADS3 * 4];


    bottom = bottomd;
    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    if (k < bottom) k += inc;

    register int restart = k;
    for (j = 0; j < 5; j++) {  // wait-free pre-passes
      // iterate over all cells assigned to thread
        while (k <= nnodesd) {
            if (massd[k] < 0) {
                for (i = 0; i < 4; i++) {
                    ch = childd[k * 4 + i];
                    child[i * THREADS3 + threadIdx.x] = ch;  // cache children 
                    if ((ch >= nbodiesd) && (massd[ch] < 0)) {
                        break;
                    }
                }

                if (i == 4) {
                    // all children are ready
                    cm = 0;
                    cg = 0;
                    cnt = 0;
#ifdef USE_DIP
                    cen.x = posd[k].x;
                    cen.y = posd[k].y;
                    cd.x = cd.y = 0;
#ifdef USE_QUA
                    cq.x = cq.y = 0;
#ifdef USE_OCT
                    co.x = co.y = 0;
#ifdef USE_HEX
                    chex.x = chex.y = 0;
#endif
#endif
#endif
#endif
                    for (i = 0; i < 4; i++) {
                        //ch = childd[k * 4 + i];
                        ch = child[i * THREADS3 + threadIdx.x];
                        if (ch >= 0) {
                            g = momd[ch].gam;
#ifdef USE_DIP
                            dr.x = posd[ch].x - cen.x; dr.y = posd[ch].y - cen.y;
                            d.x = momd[ch].dip.x; d.y = momd[ch].dip.y;
#ifdef USE_QUA

                            q.x = momd[ch].qua.x; q.y = momd[ch].qua.y;
#ifdef USE_OCT

                            o.x = momd[ch].oct.x; o.y = momd[ch].oct.y;
#ifdef USE_HEX

                            h.x = momd[ch].hex.x; h.y = momd[ch].hex.y;
#endif
#endif
#endif
#endif
                            m = massd[ch];

                            cnt += (ch >= nbodiesd) ? countd[ch] : 1;

                            // add child's contribution                                                        
                            cg += g;
#ifdef USE_DIP
                            mh = g * dr;
                            cd += d + mh;
#ifdef USE_QUA
                            mh = multz(mh, dr);
                            dh = multz(d, dr);
                            cq += q + 2 * dh + mh;
#ifdef USE_OCT
                            mh = multz(mh, dr);
                            dh = multz(dh, dr);
                            qh = multz(q, dr);
                            co += o + 3 * qh + 3 * dh + mh;
#ifdef USE_HEX
                            mh = multz(mh, dr);
                            dh = multz(dh, dr);
                            qh = multz(qh, dr);
                            oh = multz(o, dr);
                            chex += h + 4 * oh + 6 * qh + 4 * dh + mh;
#endif
#endif
#endif
#endif

                            cm += m;
                        }
                    }
                    countd[k] = cnt;

                    momd[k].gam = cg;
#ifdef USE_DIP
                    momd[k].dip.x = cd.x;
                    momd[k].dip.y = cd.y;
#ifdef USE_QUA
                    momd[k].qua.x = cq.x;
                    momd[k].qua.y = cq.y;
#ifdef USE_OCT
                    momd[k].oct.x = co.x;
                    momd[k].oct.y = co.y;
#ifdef USE_HEX
                    momd[k].hex.x = chex.x;
                    momd[k].hex.y = chex.y;
#endif
#endif
#endif
#endif

                    __threadfence();  // make sure data are visible before setting mass

                    massd[k] = cm;
                }
            }
            k += inc;  // move on to next cell
        }
        k = restart;
    }

    flag = 0;
    j = 0;
    // iterate over all cells assigned to thread
    while (k <= nnodesd) {
        if (massd[k] >= 0) {
            k += inc;
        }
        else {
            if (j == 0) {
                j = 4;
                for (i = 0; i < 4; i++) {
                    ch = childd[k * 4 + i];
                    child[i * THREADS3 + threadIdx.x] = ch;  // cache children                                 

                    if ((ch < nbodiesd) || (massd[ch] >= 0)) {
                        j--;
                    }

                }
            }
            else {
                j = 4;
                for (i = 0; i < 4; i++) {
                    //ch = childd[k * 4 + i];
                    ch = child[i * THREADS3 + threadIdx.x];

                    if ((ch < nbodiesd) || (massd[ch] >= 0))
                    {
                        j--;
                    }
                }
            }

            if (j == 0) {
                // all children are ready
                cg = 0;
                cm = 0;
#ifdef USE_DIP
                cen.x = posd[k].x;
                cen.y = posd[k].y;
                cd.x = cd.y = 0;
#ifdef USE_QUA
                cq.x = cq.y = 0;
#ifdef USE_OCT
                co.x = co.y = 0;
#ifdef USE_HEX
                chex.x = chex.y = 0;
#endif
#endif
#endif
#endif

                cnt = 0;
                for (i = 0; i < 4; i++) {
                    //ch = childd[k * 4 + i];
                    ch = child[i * THREADS3 + threadIdx.x];

                    if (ch >= 0) {
                        g = momd[ch].gam;
#ifdef USE_DIP
                        dr.x = posd[ch].x - cen.x; dr.y = posd[ch].y - cen.y;
                        d.x = momd[ch].dip.x; d.y = momd[ch].dip.y;
#ifdef USE_QUA

                        q.x = momd[ch].qua.x; q.y = momd[ch].qua.y;
#ifdef USE_OCT

                        o.x = momd[ch].oct.x; o.y = momd[ch].oct.y;
#ifdef USE_HEX

                        h.x = momd[ch].hex.x; h.y = momd[ch].hex.y;
#endif
#endif
#endif
#endif
                        m = massd[ch];

                        cnt += (ch >= nbodiesd) ? countd[ch] : 1;

                        // add child's contribution                                                        
                        cg += g;
#ifdef USE_DIP
                        mh = g * dr;
                        cd += d + mh;
#ifdef USE_QUA
                        mh = multz(mh, dr);
                        dh = multz(d, dr);
                        cq += q + 2 * dh + mh;
#ifdef USE_OCT
                        mh = multz(mh, dr);
                        dh = multz(dh, dr);
                        qh = multz(q, dr);
                        co += o + 3 * qh + 3 * dh + mh;
#ifdef USE_HEX
                        mh = multz(mh, dr);
                        dh = multz(dh, dr);
                        qh = multz(qh, dr);
                        oh = multz(o, dr);
                        chex += h + 4 * oh + 6 * qh + 4 * dh + mh;
#endif
#endif
#endif
#endif
                        cm += m;
                    }
                }
                countd[k] = cnt;

                momd[k].gam = cg;
#ifdef USE_DIP
                momd[k].dip.x = cd.x;
                momd[k].dip.y = cd.y;
#ifdef USE_QUA
                momd[k].qua.x = cq.x;
                momd[k].qua.y = cq.y;
#ifdef USE_OCT
                momd[k].oct.x = co.x;
                momd[k].oct.y = co.y;
#ifdef USE_HEX
                momd[k].hex.x = chex.x;
                momd[k].hex.y = chex.y;
#endif
#endif
#endif
#endif
                flag = 1;
            }
        }

        __threadfence();
        __syncthreads();

        if (flag != 0) {

            //atomicExch((real*)&massd[k], cm);
            massd[k] = cm;
            k += inc;
            flag = 0;
        }
    }
}





/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS4, FACTOR4)
void SortKernel(int nnodesd, int nbodiesd, volatile int* __restrict sortd, const int* __restrict countd, volatile int* __restrict startd, volatile int* __restrict childd)
{
    register int i, j, k, ch, dec, start, bottom;

    bottom = bottomd;
    dec = blockDim.x * gridDim.x;
    k = nnodesd + 1 - dec + threadIdx.x + blockIdx.x * blockDim.x;

    // iterate over all cells assigned to thread
    while (k >= bottom) {
        start = startd[k];
        if (start >= 0) {
            j = 0;
            for (i = 0; i < 4; i++) {
                ch = childd[k * 4 + i];
                if (ch >= 0) {
                    if (i != j) {
                        // move children to front (needed later for speed)
                        childd[k * 4 + i] = -1;
                        childd[k * 4 + j] = ch;
                    }
                    j++;
                    if (ch >= nbodiesd) {
                        // child is a cell
                        startd[ch] = start;  // set start ID of child
                        start += countd[ch];  // add #bodies in subtree
                    }
                    else {
                        // child is a body
                        sortd[start] = ch;  // record body in 'sorted' array
                        start++;
                    }
                }
            }
            k -= dec;  // move on to next cell
        }
    }
}


/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/


__global__
__launch_bounds__(THREADS5, FACTOR5)
void ForceCalculationKernel(int nnodesd, int nbodiesd, int* __restrict errd, real itolsqd, real epssqd, const int* __restrict sortd,
    const int* __restrict childd, const moms* __restrict momd,
    const real2* __restrict posd, volatile real2* __restrict veld)

{
    register int i, j, k, n, depth, base, sbase, diff, pd, nd;
    register real2 p, v, dr;
    register real r2;
    register moms mom;
#ifdef USE_DIP
    register real2 th;
#endif
    __shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE], node[MAXDEPTH * THREADS5 / WARPSIZE];
    __shared__ real dq[MAXDEPTH * THREADS5 / WARPSIZE];

    if (0 == threadIdx.x) {
        r2 = radiusd * 2;
        // precompute values that depend only on tree level
        dq[0] = r2 * r2;
        for (i = 1; i < maxdepthd; i++) {
            dq[i] = dq[i - 1] / 4;
            //dq[i - 1];// += epssqd;
        }
        //dq[i - 1];// += epssqd;

        if (maxdepthd > MAXDEPTH) {
            *errd = maxdepthd;
        }
    }
    __syncthreads();

    if (maxdepthd <= MAXDEPTH) {
        // figure out first thread in each warp (lane 0)
        base = threadIdx.x / WARPSIZE;
        sbase = base * WARPSIZE;
        j = base * MAXDEPTH;
        diff = threadIdx.x - sbase;
        // make multiple copies to avoid index calculations later
        if (diff < MAXDEPTH) {
            dq[diff + j] = dq[diff];
        }
        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x) {
            i = sortd[k];  // get permuted/sorted index
            // cache position info
            p.x = posd[i].x;
            p.y = posd[i].y;

            v.x = 0;
            v.y = 0;

            // initialize iteration stack, i.e., push root node onto stack
            depth = j;
            if (sbase == threadIdx.x) {
                pos[j] = 0;
                node[j] = nnodesd * 4;
            }

            do {
                // stack is not empty
                pd = pos[depth];
                nd = node[depth];
                while (pd < 4) {
                    // node on top of stack has more children to process
                    n = childd[nd + pd];  // load child pointer
                    pd++;

                    if (n >= 0) {
                        dr.x = p.x - posd[n].x;
                        dr.y = p.y - posd[n].y;
                        mom = momd[n];

                        r2 = dr.x * dr.x + dr.y * dr.y;   // compute distance squared (plus softening)
                        if ((n < nbodiesd) || __all_sync(0xffffffff, (dq[depth] + epssqd) * itolsqd < r2)) {  // check if all threads agree that cell is far enough away (or is a body)
                            real f = mom.gam / realmax(r2, epssqd);

                            v += f * dr;
#ifdef USE_DIP 
                            real2 cftr = ((r2 > 0) ? 1 / r2 : 0) * dr;

                            th = multz(cftr, cftr);
                            v += multzA(th, mom.dip);
#ifdef USE_QUA 
                            th = multz(th, cftr);
                            v += multzA(th, mom.qua);
#ifdef USE_OCT 
                            th = multz(th, cftr);
                            v += multzA(th, mom.oct);
#ifdef USE_HEX 
                            th = multz(th, cftr);
                            v += multzA(th, mom.hex);
#endif
#endif
#endif
#endif

                        }
                        else {
                            // push cell onto stack
                            if (sbase == threadIdx.x) {  // maybe don't push and inc if last child
                                pos[depth] = pd;
                                node[depth] = nd;
                            }
                            depth++;
                            pd = 0;
                            nd = n * 4;
                        }
                    }
                    else {
                        pd = 4;  // early out because all remaining children are also zero
                    }
                }
                depth--;  // done with this level
            } while (depth >= j);


            // update velocity
            veld[i].x = -v.y;
            veld[i].y = v.x;
        }
    }
}


/******************************************************************************/
/*** compute force (direct) ***************************************************/
/******************************************************************************/


__global__
//__launch_bounds__(THREADSD, FACTORD)
void ForceDirectCalculationKernel(int nnodesd, int nbodiesd,
    int* __restrict errd,
    real itolsqd, real epssqd,
    const int* __restrict sortd, const int* __restrict childd,
    const moms* __restrict momd,
    const real2* __restrict posd,
    volatile real2* __restrict veld)
{
    //*
    __shared__ real2 shr[BLOCKD];
    __shared__ real shg[BLOCKD];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    real2 pt;
    pt.x = posd[i].x;
    pt.y = posd[i].y;

    real2 vel;
    vel.x = vel.y = 0;

    real2 dr;
    real dr2, izn;


    //vortices
    for (size_t j = 0; j < nbodiesd; j += BLOCKD)
    {
        shr[threadIdx.x].x = posd[(j + threadIdx.x)].x;
        shr[threadIdx.x].y = posd[(j + threadIdx.x)].y;
        shg[threadIdx.x] = momd[(j + threadIdx.x)].gam;

        __syncthreads();

        for (size_t q = 0; q < BLOCKD; ++q)
        {
            if (j + q < nbodiesd)
            {
                dr.x = pt.x - shr[q].x;
                dr.y = pt.y - shr[q].y;
                dr2 = dr.x * dr.x + dr.y * dr.y;

                izn = shg[q] / realmax(dr2, epssqd);// / CUboundDenom(dr2, eps2); //РЎРіР»Р°Р¶РёРІР°С‚СЊ РЅР°РґРѕ!!!

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
    cudaFuncSetCacheConfig(BoundingBoxKernel, cudaFuncCachePreferShared); //1
    cudaGetLastError();  // reset error value
    
    cudaFuncSetCacheConfig(TreeBuildingKernel, cudaFuncCachePreferL1);    //2
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ClearKernel1, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ClearKernel3, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(SummarizationKernel, cudaFuncCachePreferL1);   //3
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(SortKernel, cudaFuncCachePreferL1);            //4
    cudaGetLastError();  // reset error value
    
    cudaFuncSetCacheConfig(ForceCalculationKernel, cudaFuncCachePreferEqual); //5
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ForceDirectCalculationKernel, cudaFuncCachePreferEqual); //d
    cudaGetLastError();  // reset error value
}


/******************************************************************************/





//////////////////
/// Wrappers
//////////////////



    /******************************************************************************/
    /*** initialize memory ********************************************************/
    /******************************************************************************/

    float cuInitializationKernel(int* __restrict errd)
    {
        //fprintf(stderr, "IKKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);

        cudaEventRecord(start, 0);
        InitializationKernel<<<1, 1>>> (errd);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        CudaTest("kernel 0 launch failed");
        
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/

    float cuBoundingBoxKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd, volatile int* __restrict childd,
        volatile int* __restrict massd,
        volatile moms * __restrict momd,
        volatile realPoint* __restrict posd,
        volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd)
    {
        //fprintf(stderr, "BBKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        BoundingBoxKernel<<<blocks * FACTOR1, THREADS1>>> (nnodesd, nbodiesd, startd, childd, massd, momd, (real2*)posd, (real2*)maxrd, (real2*)minrd);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 1 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/

    float cuClearKernel1(int nnodesd, int nbodiesd, volatile int* __restrict childd)
    {
        //fprintf(stderr, "C1Kernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);        

        ClearKernel1 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, childd);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel clear1 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    float cuTreeBuildingKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict errd, volatile int* __restrict childd,
        volatile realPoint* __restrict posd)
    {
        //fprintf(stderr, "TBKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        TreeBuildingKernel << <blocks * FACTOR2, THREADS2 >> > (nnodesd, nbodiesd, errd, childd, (real2*)posd);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);        
        return time;
    }

    float cuClearKernel23(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd,
        volatile int* __restrict massd,
        const real* __restrict gamd,
        volatile moms* __restrict momd)
    {
        //fprintf(stderr, "CxKernel\n");
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ClearKernel2 << <blocks * 1, 1024 >> > (nnodesd, startd, massd);
        ClearKernel3 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, gamd, momd);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel clear2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    /******************************************************************************/
    /*** compute multipole moments for all the cells ******************************/
    /******************************************************************************/

    float cuSummarizationKernel(
        const int nnodesd, const int nbodiesd,
        volatile int* __restrict countd, const int* __restrict childd,
        volatile int* __restrict massd,
        volatile moms * __restrict momd,
        volatile realPoint* __restrict posd)
    {
        //fprintf(stderr, "SKKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        SummarizationKernel << <blocks * FACTOR3, THREADS3 >> > (nnodesd, nbodiesd, countd, childd, massd, momd, (real2*)posd);
        
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        
        CudaTest("kernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }

    /******************************************************************************/
    /*** sort bodies **************************************************************/
    /******************************************************************************/

    float cuSortKernel(
        int nnodesd, int nbodiesd,
        volatile int* __restrict sortd, const int* __restrict countd,
        volatile int* __restrict startd, volatile int* __restrict childd)
    {
        //fprintf(stderr, "SRKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);        
        cudaEventRecord(start, 0);
                
        SortKernel << <blocks * FACTOR4, THREADS4 >> > (nnodesd, nbodiesd, sortd, countd, startd, childd);
        
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        
        CudaTest("kernel 4 launch failed");
 
        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


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
        volatile realPoint* __restrict veld)
    {
        //fprintf(stderr, "FCKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        
        ForceCalculationKernel << <blocks * FACTOR5, THREADS5 >> > (nnodesd, nbodiesd, errd, itolsqd, epssqd, sortd, childd, momd, (real2*)posd, (real2*)veld);
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
        int* __restrict errd,
        real itolsqd, real epssqd,
        const int* __restrict sortd, const int* __restrict childd,
        const moms * __restrict momd,
        const realPoint * __restrict posd,
        volatile realPoint* __restrict veld)
    {
        //fprintf(stderr, "DFKernel\n");
        
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);
        
        ForceDirectCalculationKernel<<<(nbodiesd + BLOCKD - 1) / BLOCKD, BLOCKD>>> (nnodesd, nbodiesd, errd, itolsqd, epssqd, sortd, childd, momd, (real2*)posd, (real2*)veld);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);
        
        CudaTest("kernel direct launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);
        
        return timeD;
    }

}//namespace BHcu


