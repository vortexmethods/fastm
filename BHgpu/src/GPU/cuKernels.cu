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
\version 1.4
\date 28 марта 2023 г.
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
#define MAXDEPTH 32
#define BLOCKD 32

#define codeLength 14
#define twoPowCodeLength (1 << codeLength)

namespace BHcu
{
    int blocks;
__device__ volatile int bottomd, maxdepthd;
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
    maxdepthd = 1;
    blkcntd = 0;
}


/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS1, FACTOR1)
void MBoundingBoxKernel(int nbodiesd, volatile real2* __restrict posd, volatile real2* __restrict Mposd, volatile real2* __restrict maxrd, volatile real2* __restrict minrd)
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
            k = 0;

            Mposd[k].x = (minr.x + maxr.x) / 2;
            Mposd[k].y = (minr.y + maxr.y) / 2;            
        }
    }
}

/******************************************************************************/
/*** Morton codes *************************************************************/
/******************************************************************************/

__global__
void MMortonCodesKernel (int nbodies, volatile real2* __restrict posd, 
    volatile int* __restrict MmortonCodesKeyUnsortd, volatile int* __restrict MmortonCodesIdxUnsortd)
{
	int bdy = blockDim.x * blockIdx.x + threadIdx.x;

	if (bdy < nbodies)
	{
		real x = twoPowCodeLength * posd[bdy].x;
		real y = twoPowCodeLength * posd[bdy].y;

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
    int nbodies, 
    int* __restrict MmortonCodesKeyd, 
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


        int gamma = i + s * d +    d * (d < 0);   //последнее слагаемое = std::min(d, 0);

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
    int nbodies,
    int* __restrict MmortonCodesKeyd,
    real2* __restrict Mposd,
    real2* __restrict Msized,
    int2* __restrict Mranged,
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
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(1024, 1)
void ClearKernel2(int nnodesd, int nbodiesd, volatile int* __restrict startd, volatile int* __restrict massd)
{
    register int k, inc, bottom;

    bottom = nnodesd - (nbodiesd - 2); //bottomd;
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
void ClearKernel3(int nnodesd, int nbodies, const real* __restrict gamd, volatile real* __restrict momsd)
{
    register int k, inc;

    inc = blockDim.x * gridDim.x;
    k = threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    
    // iterate over all cells assigned to thread
    while (k < nnodesd) {   
          momsd[k * (order * 2 - 1)] = (k < nbodies) ? gamd[k] : 0;
          for(int s = 1; s < (order * 2 - 1); ++s)
            momsd[k * (order * 2 - 1) + s] = 0;

        k += inc;
    }
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS3, FACTOR3)
void SummarizationKernel2(const int nnodesd, const int nbodiesd, volatile int* __restrict countd, const int* __restrict childd, volatile int* __restrict massd, volatile real* __restrict momsd, volatile real2* __restrict posd,
    const int* __restrict cftl)
{
    register int i, j, k, ch, inc, cnt, bottom, flag;

    register real mom[order * 2 - 1];
    //register volatile real* mom;
    
    register real cmom[order * 2 - 1];
    
    register real2 momh[order-1]; // (order-1) !!!
    register real2 cen, dr;

    register int m, cm;

    register real binCft;

    bottom = nnodesd - (nbodiesd - 2);// bottomd;
      
    __syncthreads();

    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    
    if (k < bottom) k += inc;

    //register int restart = k;

    flag = 0;
    j = 0;
    // iterate over all cells assigned to thread
    while (k <= nnodesd) {
        if (massd[k] >= 0) {
            k += inc;
        }
        else {
            if (j == 0) {
                j = 2;
                for (i = 0; i < 2; i++) {
                    ch = childd[k * 2 + i];
                    if ((ch < nbodiesd) || (massd[ch] >= 0)) {
                        j--;
                    }

                }
            }
            else {
                j = 2;
                for (i = 0; i < 2; i++) {
                    ch = childd[k * 2 + i];                    

                    if ((ch < nbodiesd) || (massd[ch] >= 0))
                    {
                        j--;
                    }
                }
            }

            if (j == 0) {
                // all children are ready

                cm = 0;

                for (int s = 0; s < (order * 2 - 1); ++s)
                    cmom[s] = 0;

                if (order > 1)
                {
                    cen.x = posd[k].x;
                    cen.y = posd[k].y;
                }


                cnt = 0;
                for (i = 0; i < 2; i++) {
                    ch = childd[k * 2 + i];                   

                    if (ch >= 0) {
                        for (int s = 0; s < (order * 2 - 1); ++s)
                            mom[s] = momsd[ch * (order * 2 - 1) + s];
                        //mom = momsd + (ch * (order * 2 - 1));

                        if (order > 1)
                        {
                            dr.x = posd[ch].x - cen.x; 
                            dr.y = posd[ch].y - cen.y;
                        }
                        
                        m = massd[ch];

                        cnt += (ch >= nbodiesd) ? countd[ch] : 1;

                        // add child's contribution    

                        cmom[0] += mom[0];

                        for (int p = 1; p < 2 * (order - 1); p += 2)
                        {
                            cmom[p + 0] += mom[p + 0];
                            cmom[p + 1] += mom[p + 1];

                            for (int q = 0; q < (p - 1) / 2; ++q)
                                momh[q] = multz(momh[q], dr);                                

                            if (p==1)
                              momh[0] = mom[0] * dr;
                            else
                              momh[(int)(p - 1) / 2] = multz(mom[p - 2], mom[p - 1], dr);
                            
                            cmom[p + 0] += momh[0].x;
                            cmom[p + 1] += momh[0].y;

                            for (int q = 1; q < (p + 1) / 2; ++q)
                            //for (int q = (p + 1) / 2 - 1; q>=0; --q)
                            {
                                binCft = cftl[((p + 1) / 2) * order + q]; 
                                cmom[p + 0] += binCft * momh[q].x;
                                cmom[p + 1] += binCft * momh[q].y;
                            }
                            
                            
                        }

                        cm += m;
                    }
                }
                countd[k] = cnt;

                for (int s = 0; s < (order * 2 - 1); ++s)
                    momsd[k * (order * 2 - 1) + s] = cmom[s];

                flag = 1;
            }
        }

        __threadfence();
        __syncthreads();

        if (flag != 0) {

            atomicExch((int*)&massd[k], cm);
            //massd[k] = cm;

            k += inc;
            flag = 0;
        }
    }

    //for (int k = 0; k <= nnodesd; ++k)
    //if ((countd[k] != massd[k]) && (massd[k] >= 0))
    
    //for (int k = 0; k <= nnodesd; ++k)
    //if (massd[k] >= 0)
    //    if (countd[k]!=0 || massd[k]!=1)
    //        if (countd[k] != massd[k])
    //            printf("k = %d, count = %d, massd = %d\n", k, countd[k], massd[k]);

    for (int k = 0; k <= nnodesd; ++k)
        if ((massd[k] < 0) && (countd[k]!=0))
            printf("k = %d, count = %d, massd = %d\n", k, countd[k], massd[k]);
}



/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS4, FACTOR4)
void SortKernel2(int nnodesd, int nbodiesd, volatile int* __restrict sortd, const int* __restrict countd, volatile int* __restrict startd, volatile int* __restrict childd)
{
    register int i, j, k, ch, dec, start, bottom;

    bottom = nnodesd - (nbodiesd - 2); //bottomd;

    //printf("startd[-1] = %d\n", startd[nnodesd]);

    dec = blockDim.x * gridDim.x;
    k = nnodesd + 1 - dec + threadIdx.x + blockIdx.x * blockDim.x;

    // iterate over all cells assigned to thread
    while (k >= bottom) {
        start = startd[k];
        if (start >= 0) {
            j = 0;
            for (i = 0; i < 2; i++) {
                ch = childd[k * 2 + i];
                if (ch >= 0) {
                    if (i != j) {
                        // move children to front (needed later for speed)
                        //printf("sort: k = %d, childd = %d, %d, ch = %d, i = %d, j= %d\n", k, childd[k * 2 + 0], childd[k * 2 + 1], ch, i, j);
                        childd[k * 2 + i] = -1;
                        childd[k * 2 + j] = ch;
                        
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
void ForceCalculationKernel2(int nnodesd, int nbodiesd, \
    int* __restrict errd, real itolsqd, real epssqd, const int* __restrict sortd,
    const int* __restrict childd, const real* __restrict momsd,
    const real2* __restrict posd, volatile real2* __restrict veld, 
    const real2* __restrict Msized)

{
register int i, j, k, n, depth, base, sbase, /*diff,*/ pd, nd;
register real2 p, v, dr;
register real r2;
//register real mom[(order * 2 - 1)];
register const real* mom;

register real2 th;

__shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE], node[MAXDEPTH * THREADS5 / WARPSIZE];

maxdepthd = 28; ////////!!!!!!!!

if (maxdepthd <= MAXDEPTH)
{
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
        i = sortd[k];  // get permuted/sorted index
        //if (i != k)
        //    printf("i = %d, k = %d\n", i, k);


        // cache position info
        p.x = posd[i].x;
        p.y = posd[i].y;

        v.x = 0;
        v.y = 0;

        // initialize iteration stack, i.e., push root node onto stack
        depth = j;
        if (sbase == threadIdx.x)
        {
            pos[j] = 0;
            node[j] = nnodesd * 2;
        }

        do
        {
            // stack is not empty
            pd = pos[depth];
            nd = node[depth];

            while (pd < 2)
            {
                // node on top of stack has more children to process
                n = childd[nd + pd];  // load child pointer
                ++pd;

                if (n >= 0)
                {
                    dr = p - posd[n];

                    mom = momsd + (n * (order * 2 - 1));

                    r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared (plus softening)

                    real sumSide = Msized[n].x + Msized[n].y;

                    if ((n < nbodiesd) || __all_sync(0xffffffff, (sumSide * sumSide + epssqd) * itolsqd < r2))
                        {  // check if all threads agree that cell is far enough away (or is a body)

                            real f = mom[0] / realmax(r2, epssqd);
                            v += f * dr;

                            if (order > 1)
                            {
                                real2 cftr = (r2 ? (1 / r2) : 0) * dr;

                                //s=1:
                                th = multz(cftr, cftr);
                                v += multzA(th, mom[1], mom[2]);

                                for (int s = 3; s < (order * 2 - 1); s += 2)
                                {
                                    th = multz(th, cftr);
#ifdef CALCinFLOAT                                    
                                    if (isinf(th.x) || isinf(th.y))
                                    {
                                        //printf("s = %d\n", s);
                                        break;
                                    }
#endif
                                    v += multzA(th, mom[s], mom[s + 1]);
                                }
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
                            nd = n * 2;
                        }
                    }
                    else
                    {
                        pd = 2;  // early out because all remaining children are also zero
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
    const real * __restrict momsd,
    const real2* __restrict posd,
    volatile real2* __restrict veld)
{
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
        shg[threadIdx.x] = momsd[(j + threadIdx.x)*(order * 2 -1)];//momd[(j + threadIdx.x)].gam;

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

    cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ClearKernel3, cudaFuncCachePreferL1);
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
	float McuBoundingBoxKernel(
		int nbodiesd,
		volatile realPoint* __restrict posd,
		volatile realPoint* __restrict Mposd,
		volatile realPoint* __restrict maxrd, volatile realPoint* __restrict minrd)
	{
		cudaEvent_t start, stop;
		float time;

		cudaEventCreate(&start);  cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		MBoundingBoxKernel<<<blocks * FACTOR1, THREADS1>>> (nbodiesd, (real2*)posd, (real2*)Mposd, (real2*)maxrd, (real2*)minrd);
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
        realPoint* __restrict posd,
        int* __restrict MmortonCodesKeyUnsortd, int* __restrict MmortonCodesIdxUnsortd,
        int* __restrict MmortonCodesKeyd, int* __restrict MmortonCodesIdxd,
        intPair* __restrict Mranged
        )
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = (nbodiesd + 31) / 32;
        dim3 Mthreads = 32;

        MMortonCodesKernel << <Mblocks, Mthreads >> > (nbodiesd, (real2*)posd, MmortonCodesKeyUnsortd, MmortonCodesIdxUnsortd);


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
        int* __restrict MmortonCodesKeyd,
        realPoint* __restrict Mposd,
        realPoint* __restrict Msized,
        intPair* __restrict Mranged,
        int* __restrict MlevelUnsortd,
        int* __restrict MlevelSortd,
        int* __restrict MindexUnsortd,
        int* __restrict MindexSortd
    )
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((nbodiesd - 1) + 31) / 32;
        dim3 Mthreads = 32;

		MMortonInternalCellsGeometryKernel << <Mblocks, Mthreads >> > (nbodiesd, MmortonCodesKeyd, (real2*)Mposd, (real2*)Msized, (int2*)Mranged,
            MlevelUnsortd, MindexUnsortd);


        RadixSortFromCUB( \
            MlevelUnsortd, MlevelSortd, \
            MindexUnsortd, MindexSortd, \
            nbodiesd-1, 0, 2 * codeLength);


        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }



	 	 

    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/   

    float cuClearKernel23(
        int nnodesd, int nbodiesd,
        volatile int* __restrict startd,
        volatile int* __restrict massd,
        const real* __restrict gamd,
        volatile real* __restrict momsd)
    {
        //fprintf(stderr, "CxKernel\n");
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ClearKernel2 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, startd, massd);
        ClearKernel3 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, gamd, momsd);

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
        volatile int* __restrict countd, const int* __restrict childd,
        volatile int* __restrict massd,
        volatile real* __restrict momsd,
        volatile realPoint* __restrict posd,
        const int* __restrict cftl)
    {
        //fprintf(stderr, "SKKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        SummarizationKernel2 << <blocks * FACTOR3, THREADS3 >> > (nnodesd, nbodiesd, countd, childd, massd, momsd, (real2*)posd, cftl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }

    /******************************************************************************/
    /*** sort bodies **************************************************************/
    /******************************************************************************/
    float cuSortKernel2(
        int nnodesd, int nbodiesd,
        volatile int* __restrict sortd, const int* __restrict countd,
        volatile int* __restrict startd, volatile int* __restrict childd)
    {
        //fprintf(stderr, "SRKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        SortKernel2 << <blocks * FACTOR4, THREADS4 >> > (nnodesd, nbodiesd, sortd, countd, startd, childd);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 4 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }

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
        volatile realPoint* __restrict Msized)
    {
        //fprintf(stderr, "FCKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ForceCalculationKernel2 << <blocks * FACTOR5, THREADS5 >> > (nnodesd, nbodiesd, errd, itolsqd, epssqd, sortd, childd, momsd, (real2*)posd, (real2*)veld, (real2*)Msized);
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
        const real * __restrict momsd,
        const realPoint * __restrict posd,
        volatile realPoint* __restrict veld)
    {
        //fprintf(stderr, "DFKernel\n");
        
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);
        
        ForceDirectCalculationKernel<<<(nbodiesd + BLOCKD - 1) / BLOCKD, BLOCKD>>> (nnodesd, nbodiesd, errd, itolsqd, epssqd, sortd, childd, momsd, (real2*)posd, (real2*)veld);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);
        
        CudaTest("kernel direct launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);
        
        return timeD;
    }

}//namespace BHcu