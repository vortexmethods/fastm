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
| File name: main.cpp                                                         |
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
 \brief Barnes-Hut method (CUDA) for 2D vortex particles + Morton tree
 \author Марчевский Илья Константинович
 \author Рятина Евгения Павловна
 \author Колганова Александра Олеговна
 \version 1.4
 \date 28 марта 2023 г.
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

//Контроль правильности построения дерева
    void traverseT(int v, int* TEMPchild)
    {
        int v0 = TEMPchild[2*v+0];
        int v1 = TEMPchild[2*v+1];
        if (v0 > 0)
            traverseT(v0, TEMPchild);
        if (v1 > 0)
            traverseT(v1, TEMPchild);    
    } 


int main(int argc, char** argv)
{
	/// ПРОКОММЕНТИРОВАННАЯ ЧАСТЬ
	//без буквы "l" на конце - на host,
	//c буквой "l" на конце - на device

	std::vector<real> gam;          //циркуляции из файла
	real* gaml;
	
	std::vector<realPoint> pos;     //положения из файла
	realPoint* posl;
	
	std::vector<int> mass;          //массы (единица для точечного вихря, число вихрей для ячейки) //todo зачем они на хосте?
	int* massl;
	
	
	std::vector<realPoint> vel;     //для вычисляемых скоростей
	realPoint* vell;                // - быстрым методом
	realPoint* vellBS;              // - прямым методом Био--Савара

	std::vector<int> cft;           //биномиальные коэффициенты
	int* cftl;

	realPoint* maxrl, *minrl;       //габаритный прямоугольник
	
	real* momsl;                    //мультипольные моменты всех ячеек; хранятся в виде <mom_0, mom_1x, mom_1y, ..., mom_px, mom_py>, <для второй ячейки> ...

	int* childl;                    //номера ячеек-потомков (хранятся парами)

	//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
	int blocks;

	//техническая
	std::vector<int> error(1);
	int* errl;

	//Число ячеек дерева и тел
	int nnodes, nbodies;

	//Радиус вихря и параметр близости и их квадраты
	real epssq = (real)(EPS * EPS);
	real itolsq = (real)(1 / (THETA * THETA));

	//Статистика
	float timing[7], mintiming[7] = { 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[7]{};
	double runtime, minruntime = 1e+10, avruntime = 0;


	/// НИЖЕ пока еще не прокомментирвоанная часть
	
	
	PrintLogoToStream(std::cout);
		
	//Проверка функционирования видеокарты
	CudaSelect(dev);
	setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)
	

	   
	
    
	//Массивы исходных данных и результата на host


    //Указатели на массивы, хранящиеся на device
	int * sortl, * countl, * startl;
	
	//For Morton tree
	int* MmortonCodesKeyUnsortl;
	int* MmortonCodesIdxUnsortl;	
	int* MmortonCodesKeyl;
	int* MmortonCodesIdxl;

	int* MlevelUnsortl;
	int* MlevelSortl;
	int* MindexUnsortl;
	int* MindexSortl;

	realPoint* Mposl;
	realPoint* Msizel;
	int* Mparentl;
	intPair* Mchildl;
	intPair* Mrangel;

	int* Mflagl;
	int* Mmassl;	
	int* Mlockl;
	real* Mmoml;
	real* Mgaml;

	


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

			momsl = (real*)cudaNew((nnodes + 1) * (order * 2 - 1), sizeof(real));
			posl = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
			vell = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
			vellBS = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));

			countl = (int*)cudaNew(nnodes + 1, sizeof(int));
			startl = (int*)cudaNew(nnodes + 1, sizeof(int));
			sortl = (int*)cudaNew(nnodes + 1, sizeof(int));

			maxrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
			minrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));

			///For MortonTree
			MmortonCodesKeyUnsortl = (int*)cudaNew(nnodes + 1, sizeof(int));
			MmortonCodesIdxUnsortl = (int*)cudaNew(nnodes + 1, sizeof(int));
			MmortonCodesKeyl = (int*)cudaNew(nnodes + 1, sizeof(int));
			MmortonCodesIdxl = (int*)cudaNew(nnodes + 1, sizeof(int));

			Mposl = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
			Msizel = (realPoint*)cudaNew(nnodes + 1, sizeof(realPoint));
			Mparentl = (int*)cudaNew(nnodes + 1, sizeof(int));
			Mchildl = (intPair*)cudaNew(nnodes + 1, sizeof(intPair));
			Mrangel = (intPair*)cudaNew(nnodes + 1, sizeof(intPair));

			Mflagl = (int*)cudaNew(nnodes + 1, sizeof(int));
			Mmassl = (int*)cudaNew(nnodes + 1, sizeof(int));
			Mlockl = (int*)cudaNew(nnodes + 1, 2 * sizeof(int));
			Mmoml = (real*)cudaNew((nnodes + 1), (2 * order - 1) * sizeof(real));
			Mgaml = (real*)cudaNew((nnodes + 1), sizeof(real));

			MlevelUnsortl = (int*)cudaNew(nnodes + 1, sizeof(int));
			MlevelSortl = (int*)cudaNew(nnodes + 1, sizeof(int));;
			MindexUnsortl = (int*)cudaNew(nnodes + 1, sizeof(int));;
			MindexSortl = (int*)cudaNew(nnodes + 1, sizeof(int));;

			//MMchildl = (int*)cudaNew((nnodes + 1) * 2, sizeof(int));

		}

		if (run == 0)
		{
			cudaCopyVecToDevice(mass.data(), massl, nbodies, sizeof(int));
			cudaCopyVecToDevice(gam.data(), gaml, nbodies, sizeof(real));
			cudaCopyVecToDevice(pos.data(), posl, nbodies, sizeof(realPoint));
			cudaCopyVecToDevice(vel.data(), vell, nbodies, sizeof(realPoint));
			cudaCopyVecToDevice(vel.data(), vellBS, nbodies, sizeof(realPoint));
			cudaCopyVecToDevice(cft.data(), cftl, (order * order), sizeof(int));
		}

		// run timesteps (launch GPU kernels)

		double starttime, endtime;
		starttime = omp_get_wtime();


		timing[0] += cuInitializationKernel(errl);
	
		timing[1] += McuBoundingBoxKernel(nbodies, posl, Mposl, maxrl, minrl);

		//realPoint maxrh, minrh;
		//cudaCopyVecFromDevice(maxrl, &maxrh, 2, sizeof(real));
		//cudaCopyVecFromDevice(minrl, &minrh, 2, sizeof(real));
		//
		//std::cout << "x : " << minrh[0] << "..." << maxrh[0] << "\n";
		//std::cout << "y : " << minrh[1] << "..." << maxrh[1] << "\n";

		auto ttt1 = McuMortonCodesKernel(nbodies, posl, MmortonCodesKeyUnsortl, MmortonCodesIdxUnsortl, MmortonCodesKeyl, MmortonCodesIdxl, Mrangel);
		timing[2] += ttt1;

		//std::cout << "McuMortonCodesKernel = " << ttt1 << std::endl;

		//int mcd[2];
		//cudaCopyVecFromDevice(MmortonCodesKeyUnsortl + 0, mcd, 1, sizeof(int));
		//cudaCopyVecFromDevice(MmortonCodesKeyl + 0, mcd+1, 1, sizeof(int));
		//std::cout << mcd[0] << " " << mcd[1] << std::endl;		

		auto ttt2 = McuMortonInternalNodesKernel(nbodies, MmortonCodesKeyl, Mparentl, Mchildl, Mrangel);
		timing[2] += ttt2;

		//std::cout << "McuMortonInternalNodesKernel = " << ttt2 << std::endl;

		//std::vector<int> range(2*nbodies);
		//cudaCopyVecFromDevice(Mrangel, range.data(), 2*nbodies, sizeof(int));
		//std::ofstream parentFile("range.txt");
		//for (int i = 0; i < nbodies; ++i)
		//	parentFile << i << " " << range[2*i] << " " << range[2 * i + 1] << "\n";
		//parentFile.close();

		auto ttt3 = McuMortonInternalCellsGeometryKernel(nbodies, MmortonCodesKeyl, Mposl, Msizel, Mrangel,
			MlevelUnsortl, MlevelSortl, MindexUnsortl, MindexSortl);
		timing[2] += ttt3;

		//std::cout << "McuMortonInternalCellsGeometryKernel = " << ttt3 << std::endl;

		//std::vector<real> positions(2*nbodies);
		//std::vector<real> sizes(2*nbodies);
		//cudaCopyVecFromDevice(Mposl, positions.data(), 2*nbodies, sizeof(real));
		//cudaCopyVecFromDevice(Msizel, sizes.data(), 2*nbodies, sizeof(real));
		//std::ofstream parentFile2("pos.txt");
		//for (int i = 0; i < nbodies; ++i)
		//	parentFile2 << i << " " << positions[2*i] << " " << positions[2 * i + 1] << " " << sizes[2 * i] << " " << sizes[2 * i + 1] << "\n";
		//parentFile2.close();


		std::vector<int> MmortonCodesIdxSortedh(nnodes + 1);
		cudaCopyVecFromDevice(MmortonCodesIdxl, MmortonCodesIdxSortedh.data(), nbodies, sizeof(int));

		std::vector<real> TEMPgamh(nnodes+1);

		for (int i = 0; i < nbodies; ++i)
			TEMPgamh[i] = gam[MmortonCodesIdxSortedh[i]];
		real* TEMPgaml;
		TEMPgaml = (real*)cudaNew(nnodes + 1, sizeof(real));
		cudaCopyVecToDevice(TEMPgamh.data(), TEMPgaml, nbodies, sizeof(real));
				
		timing[2] += cuClearKernel23(nnodes, nbodies, startl, massl, TEMPgaml, momsl);

		// Later - to GPU
		std::vector<int> MlevelSorth(nbodies - 1), MindexSorth(nbodies - 1);
		std::vector<int> MlevelUnsorth(nbodies - 1), MindexUnsorth(nbodies - 1);

		cudaCopyVecFromDevice(MlevelSortl, MlevelSorth.data(), nbodies - 1, sizeof(int));
		cudaCopyVecFromDevice(MindexSortl, MindexSorth.data(), nbodies - 1, sizeof(int));

		cudaCopyVecFromDevice(MlevelUnsortl, MlevelUnsorth.data(), nbodies - 1, sizeof(int));
		cudaCopyVecFromDevice(MindexUnsortl, MindexUnsorth.data(), nbodies - 1, sizeof(int));

		std::vector<int> Mchildh((nbodies + 1) * 2);
		cudaCopyVecFromDevice(Mchildl, Mchildh.data(), (nbodies + 1) * 2, sizeof(int));

		std::vector<realPoint> positions(nbodies - 1);
		cudaCopyVecFromDevice(Mposl, positions.data(), 2 * (nbodies - 1), sizeof(real));

		//childl  ->  (nnodes+1)*4;
		//posl    ->  (nnodes+1);

		std::vector<int> MindexSorthT(nbodies - 1);
		for (int i = 0; i < MindexSorth.size(); ++i)
			MindexSorthT[MindexSorth[i]] = i;

		std::vector<int> TEMPchildh((nnodes + 1) * 2, -1);
		std::vector<realPoint> TEMPposh((nnodes + 1));
		std::vector<int> TMPmass((nbodies + 1), 1);

		/*
		for (int i = 0; i < nbodies - 1; ++i)
			traverse(i, Mchildh.data(), nbodies);

		std::cout << "Tree is correct!" << std::endl;
		*/		
		
		for (int i = 0; i < nbodies; ++i)
		{
			TEMPposh[i] = pos[MmortonCodesIdxSortedh[i]];
			/*
			{
				if (MmortonCodesIdxSortedh[i] == 0)
					std::cout << "proobraz(0) = " << i << std::endl;
				
				if (MmortonCodesIdxSortedh[i] == 1)
					std::cout << "proobraz(1) = " << i << std::endl;
			}
			*/
		}

		for (int i = 0; i < nbodies - 1; ++i)
		{
			TEMPposh[nnodes - i] = positions[MindexSorth[i]];

			if (Mchildh[2 * MindexSorth[i] + 0] >= nbodies)
				TEMPchildh[2 * (nnodes - i) + 0] = Mchildh[2 * MindexSorth[i] + 0] - nbodies;
			else
				TEMPchildh[2 * (nnodes - i) + 0] = (nnodes + 1) - 1 - MindexSorthT[Mchildh[2 * MindexSorth[i] + 0]];

			if (Mchildh[2 * MindexSorth[i] + 1] >= nbodies)
				TEMPchildh[2 * (nnodes - i) + 1] = Mchildh[2 * MindexSorth[i] + 1] - nbodies;
			else
				TEMPchildh[2 * (nnodes - i) + 1] = (nnodes + 1) - 1 - MindexSorthT[Mchildh[2 * MindexSorth[i] + 1]];
		}


		std::vector<realPoint> sizes(nbodies - 1, {0, 0});
		cudaCopyVecFromDevice(Msizel, sizes.data(), 2 * (nbodies - 1), sizeof(real));

		std::vector<realPoint> TEMPsizesh(nnodes+1, {0, 0});

		for (int i = 0; i < nbodies - 1; ++i)
		{
			TEMPsizesh[nnodes - i] = sizes[MindexSorth[i]];
		}


		/*
		for (int i = nnodes; i > nnodes - (nbodies - 1); --i)
			traverseT(i, TEMPchildh.data());

		std::cout << "TreeT is correct!" << std::endl;
		*/

		int* TEMPchildl;
		TEMPchildl = (int*)cudaNew(2 * (nnodes + 1), sizeof(int));
		cudaCopyVecToDevice(TEMPchildh.data(), TEMPchildl, (nnodes + 1) * 2, sizeof(int));

		realPoint* TEMPposl;
		TEMPposl = (realPoint*)cudaNew(2 * (nnodes + 1), sizeof(realPoint));
		cudaCopyVecToDevice(TEMPposh.data(), TEMPposl, (nnodes + 1) * 2, sizeof(real));

		realPoint* TEMPsizel;
		TEMPsizel = (realPoint*)cudaNew(2 * (nnodes + 1), sizeof(realPoint));
		cudaCopyVecToDevice(TEMPsizesh.data(), TEMPsizel, (nnodes + 1) * 2, sizeof(real));


		/// то, что выше ==> на видеокарту


		timing[3] += cuSummarizationKernel2(nnodes, nbodies, countl, TEMPchildl, massl, momsl, TEMPposl, cftl);


		timing[4] += cuSortKernel2(nnodes, nbodies, sortl, countl, startl, TEMPchildl);
		

		timing[5] += cuForceCalculationKernel2(nnodes, nbodies, errl, itolsq, epssq, sortl, TEMPchildl, momsl, TEMPposl, vell, TEMPsizel);
		
		
		std::vector<realPoint> TEMPvelh(nbodies);
		
		cudaCopyVecFromDevice(vell, TEMPvelh.data(), nbodies, sizeof(realPoint));

	

		for (int i = 0; i < nbodies; ++i)		
			vel[MmortonCodesIdxSortedh[i]] = TEMPvelh[i] * IDPI;
       
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
			float tBS = cuForceDirectCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, momsl, posl, vellBS); //cuForceDirectCalculationKernel(nnodes, nbodies, errl, itolsq, epssq, sortl, childl, moml, posl, veld);
			std::cout << "done!" << std::endl;
			std::cout << "Time (Biot-Savart): " << tBS << " ms" << std::endl;

			cudaCopyVecFromDevice(vellBS, veloBS.data(), nbodies, sizeof(realPoint));
		
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
	cudaDelete(vellBS);
    cudaDelete(momsl);
   
	cudaDelete(countl);
	cudaDelete(startl);
	cudaDelete(sortl);

	cudaDelete(maxrl);
	cudaDelete(minrl);

	///For Morton tree
	cudaDelete(MmortonCodesKeyUnsortl);
	cudaDelete(MmortonCodesIdxUnsortl);	
	cudaDelete(MmortonCodesKeyl);
	cudaDelete(MmortonCodesIdxl);

	cudaDelete(Mposl);
	cudaDelete(Msizel);
	cudaDelete(Mparentl);
	cudaDelete(Mchildl);
	cudaDelete(Mrangel);

	cudaDelete(Mflagl);
	cudaDelete(Mmassl);
	cudaDelete(Mlockl);
	cudaDelete(Mmoml);
	cudaDelete(Mgaml);

	cudaDelete(MlevelUnsortl);
	cudaDelete(MlevelSortl);
	cudaDelete(MindexUnsortl);
	cudaDelete(MindexSortl);
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