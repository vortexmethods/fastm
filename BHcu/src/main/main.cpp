/*--------------------------------*- BHcu -*-----------------*---------------*\
| #####   ##  ##                |                            | Version 1.1    |
| ##  ##  ##  ##   ####  ##  ## |  BHcu: Barnes-Hut method   | 2022/08/28     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: main.cu                                                          |
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

/*
 * Portions of this program were originally released under the following license
 *
 * CUDA BarnesHut v3.1: Simulation of the gravitational forces
 * in a galactic cluster using the Barnes-Hut n-body algorithm
 *
 * Copyright (c) 2013, Texas State University-San Marcos. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted for academic, research, experimental, or personal use provided that
 * the following conditions are met:
 *
 *    * Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *    * Neither the name of Texas State University-San Marcos nor the names of its
 *      contributors may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 *
 * For all other uses, please contact the Office for Commercialization and Industry
 * Relations at Texas State University-San Marcos <http://www.txstate.edu/ocir/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Martin Burtscher <burtscher@txstate.edu>
 *
 */

 /*!
 \file
 \brief Barnes-Hut method (CUDA) for 2D vortex particles
 \author Марчевский Илья Константинович
 \author Рятина Евгения Павловна
 \author Колганова Александра Олеговна
 \version 1.1
 \date 28 августа 2022 г.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "omp.h"

//#define NUMVECTORwithoutARRAY

#include "cuKernels.cuh"
#include "Logo.h"
#include "Params.h"
#include "Point2D.h"
#include "Vortex2D.h"


using namespace BHcu;

const real IDPI = (real)0.15915494309189534;

/******************************************************************************/
/******************************************************************************/
void setBinom(std::vector<int>& cft);

int main(int argc, char** argv)
{
	PrintLogoToStream(std::cout);

	//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
	int blocks;
	
	//Проверка функционирования видеокарты
	//const int dev = atoi(argv[3]);
	CudaSelect(dev);
	setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)
		
	//техническая
    std::vector<int> error(1);

	//Число ячеек дерева и тел
    int nnodes, nbodies;
       
	//Статистика
	float timing[7], mintiming[7] = { 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[7]{};
	double runtime, minruntime = 1e+10, avruntime = 0;
    
	//Массивы исходных данных и результата на host
    std::vector<int> mass;
	std::vector<real> gam;
    std::vector<realPoint> pos, vel;
	std::vector<int> cft;

    //Указатели на массивы, хранящиеся на device
	int* errl, * sortl, * childl, * countl, * startl, *cftl;
	int * massl;
	real * gaml;
	real * moms;

	realPoint* posl;
	realPoint* vell, * veld;
	realPoint* maxrl;
	realPoint* minrl;
        
	real epssq = (real)(EPS * EPS);
	real itolsq = (real)(1.0 / (THETA * THETA));	

	// generate input
	std::ifstream infile;
	infile.open(nameFile);
	
	infile >> nbodies;
	nnodes = nbodies * 2;
	if (nnodes < 1024 * blocks)
		nnodes = 1024 * blocks;
	while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
		nnodes++;
	nnodes--;

	try {
		mass.resize(nbodies);
		gam.resize(nbodies);
		pos.resize(nbodies);
		vel.resize(nbodies);

		cft.resize(order * (order + 1));
		setBinom(cft);
	}
	catch (...)
	{
		fprintf(stderr, "cannot allocate host memory\n");  exit(-1);
	}

	std::vector<realVortex> vtx(nbodies);

	for (int i = 0; i < nbodies; i++) {

		infile >> vtx[i].r()[0] >> vtx[i].r()[1] >> vtx[i].g();

		pos[i] = vtx[i].r();		

		mass[i] = 1;

		gam[i] = vtx[i].g();
		vel[i].toZero();	
	}

    //////////////////////////// 
    /// RUNS
    ////////////////////////////    

	KernelsOptimization();
	
	PrintConfiguration(nbodies);

    for (int run = 0; run < runs; run++) 
	{
        for (int i = 0; i < 6; i++) 
			timing[i] = 0;       

        // allocate memory
        if (run == 0) {				
			cftl = (int*)cudaNew(order * (order + 1), sizeof(int));                 
            errl = (int*)cudaNew(1, sizeof(int));
            childl = (int*)cudaNew((nnodes + 1) * 4, sizeof(int));
            massl = (int*)cudaNew(nnodes + 1, sizeof(int));
            gaml = (real*)cudaNew(nnodes + 1, sizeof(real));
			moms = (real*) cudaNew((nnodes + 1) * (order * 2 - 1), sizeof(real));
            posl = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
            vell = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
            veld = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
            
            countl = (int*)cudaNew(nnodes + 1, sizeof(int));
            startl = (int*)cudaNew(nnodes + 1, sizeof(int));
            sortl = (int*)cudaNew(nnodes + 1, sizeof(int));

            maxrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
            minrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
        }
        
        if (run == 0)
        {   
			cudaCopyVecToDevice(mass.data(), massl, nbodies, sizeof(int));
            cudaCopyVecToDevice(gam.data(), gaml, nbodies, sizeof(real));
            cudaCopyVecToDevice(pos.data(), posl, nbodies, sizeof(realPoint));
            cudaCopyVecToDevice(vel.data(), vell, nbodies, sizeof(realPoint));
            cudaCopyVecToDevice(vel.data(), veld, nbodies, sizeof(realPoint));
			cudaCopyVecToDevice(cft.data(), cftl, (order * order), sizeof(int));
        }
        
        // run timesteps (launch GPU kernels)

        double starttime, endtime;
        starttime = omp_get_wtime();

        
        timing[0] += cuInitializationKernel(errl);               
        
		timing[1] += cuBoundingBoxKernel(nnodes, nbodies, startl, childl, massl, moms, posl, maxrl, minrl); /*cuBoundingBoxKernel(nnodes, nbodies, startl, childl, massl, moml, posl, maxrl, minrl);*/
                
        timing[2] += cuClearKernel1(nnodes, nbodies, childl);
        timing[2] += cuTreeBuildingKernel(nnodes, nbodies, errl, childl, posl);
        timing[2] += cuClearKernel23(nnodes, nbodies, startl, massl, gaml, moms); /*cuClearKernel23(nnodes, nbodies, startl, massl, gaml, moml);*/
                
        timing[3] += cuSummarizationKernel(nnodes, nbodies, countl, childl, massl, moms, posl, cftl);/*cuSummarizationKernel(nnodes, nbodies, countl, childl, massl, moml, posl);*/
	    
        timing[4] += cuSortKernel(nnodes, nbodies, sortl, countl, startl, childl);

        timing[5] += cuForceCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, moms, posl, vell); /*cuForceCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, moml, posl, vell);*/
            
        // transfer result back to CPU
		//cudaCopyVecFromDevice(errl, error.data(), 1, sizeof(int));
		//cudaCopyVecFromDevice(posl, pos.data(), nbodies, sizeof(real2));
		cudaCopyVecFromDevice(vell, vel.data(), nbodies, sizeof(realPoint));

		for (int i = 0; i < nbodies; ++i)		
			vel[i] *= IDPI;					
       
        endtime = omp_get_wtime();	
		runtime = endtime - starttime;
		if (minruntime > runtime)
			minruntime = runtime;
		avruntime += runtime;

		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		for (int i = 1; i <= 6; ++i)
		{
			if (mintiming[i] > timing[i])
				mintiming[i] = timing[i];
			avtiming[i] += timing[i];
		}

		PrintStatistics(run, runs, error[0], timing, mintiming, avtiming, runtime, minruntime, avruntime);
    }

	if (compare)
	{
		PrintAccuracyHead();		

		std::vector<realPoint> veloBS(nbodies);

		bool realBSfromFile = BSfromFile;
		if (BSfromFile)
		{
			std::ifstream f("../../res/velBS" + task + ".txt");
			if (realBSfromFile = f.good())
				std::cout << "File with Biot-Savart results is found, loading it... " << std::flush;
			else
				std::cout << "File with Biot-Savart results is not found, computing it... " << std::flush;
		}

		if (!realBSfromFile)
		{
			float tBS = cuForceDirectCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, moms, posl, veld); //cuForceDirectCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, moml, posl, veld);
			std::cout << "done!" << std::endl;
			std::cout << "Time (Biot-Savart): " << tBS << " ms" << std::endl;

			cudaCopyVecFromDevice(veld, veloBS.data(), nbodies, sizeof(realPoint));
		
			for (int i = 0; i < nbodies; ++i)
				veloBS[i] *= IDPI;

			std::ofstream velFileBS("../../res/velBS" + task + ".txt");
			velFileBS.precision(16);
			for (int i = 0; i < nbodies; i++)
				velFileBS << veloBS[i][0] << " " << veloBS[i][1] << std::endl;
			velFileBS.close();
		}
		else
		{
			std::ifstream velFileBS("../../res/velBS" + task + ".txt");
			for (int i = 0; i < nbodies; i++) 
				velFileBS >> veloBS[i][0] >> veloBS[i][1];			
			velFileBS.close();
			std::cout << "done!" << std::endl;
		}

		real err = 0, absVel = 0;
		for (int i = 0; i < nbodies; i++) 
		{
			err += (vel[i] - veloBS[i]).length<real>();
			absVel += veloBS[i].length<real>();
		}
		
		PrintAccuracyError(err / absVel);

		//printf("BHcu: %.6e %.6e\n", vel[0][0], vel[0][1]);
		//printf("BScu: %.6e %.6e\n", veloBS[0][0], veloBS[0][1]);
	}

	if (save)
	{
		std::ofstream velFile("../../res/BH" + task + ".txt");
		velFile.precision(16);
		for (int i = 0; i < nbodies; i++)
			velFile << vel[i][0] << " " << vel[i][1] << std::endl;
		velFile.close();
	}
	std::cout << "Goodbye! " << std::endl;
    
    cudaDelete(errl);
	cudaDelete(childl);

	cudaDelete(massl);
	cudaDelete(posl);
	cudaDelete(vell);
	cudaDelete(veld);
   
	cudaDelete(countl);
	cudaDelete(startl);
	cudaDelete(sortl);

	cudaDelete(maxrl);
	cudaDelete(minrl);

    return 0;
}

void setBinom(std::vector<int>& cft)
{    
    cft[0] = 1;

    for (int i = 1; i < order; ++i)
    {
        cft[i * order + 0] = 1;
        cft[i * order + i] = 1;
        for (int j = 1; j < i; ++j)
            cft[i * order + j] = cft[(i - 1) * order + j] + cft[(i - 1) * order + (j - 1)];
    }
}