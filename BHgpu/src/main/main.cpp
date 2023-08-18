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
 \version 1.5
 \date 29 августа 2023 г.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "omp.h"

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
/*
    void traverseT(int v, int* TEMPchild)
    {
        int v0 = TEMPchild[2*v+0];
        int v1 = TEMPchild[2*v+1];
        if (v0 > 0)
            traverseT(v0, TEMPchild);
        if (v1 > 0)
            traverseT(v1, TEMPchild);    
    } 

	
	//for (int i = 0; i < nbodies - 1; ++i)
	//	traverse(i, Mchildh.data(), nbodies);
	//std::cout << "Tree is correct!" << std::endl;
*/


int main(int argc, char** argv)
{
	/// ПРОКОММЕНТИРОВАННАЯ ЧАСТЬ
	//без буквы "l" на конце - на host,
	//c буквой "l" на конце - на device	
	std::vector<realVortex> vtx;    //вихри из файла
	realVortex* vtxl;

	int* massl;	//массы (единица для точечного вихря, число вихрей для ячейки) //todo зачем они на хосте?	
	
	std::vector<realPoint> vel;     //для вычисляемых скоростей
	realPoint* vell;                // - быстрым методом
	realPoint* vellBS;              // - прямым методом Био--Савара

	std::vector<int> cft;           //биномиальные коэффициенты
	
	realPoint* maxrl, *minrl;       //габаритный прямоугольник
	
	realPoint* momsl;               //мультипольные моменты всех ячеек; хранятся в виде <mom_0x, mom_0y=0 mom_1x, mom_1y, ..., mom_px, mom_py>, <для второй ячейки> ...

	//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
	int blocks;

	//Число ячеек дерева и тел
	int nnodes, nbodies;

	//Радиус вихря и параметр близости и их квадраты
	real epssq = (real)(EPS * EPS);
	real itolsq = (real)(1 / (THETA * THETA));

	//Статистика
	float timing[7], mintiming[7] = { 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[7]{};
	double runtime, minruntime = 1e+10, avruntime = 0;

	
	PrintLogoToStream(std::cout);
		
	//Проверка функционирования видеокарты
	CudaSelect(dev);
	setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		
    
	//Массивы исходных данных и результата на host
	//Указатели на массивы, хранящиеся на device

	//For Morton tree
	int* MmortonCodesKeyUnsortl;
	int* MmortonCodesKeyl;

	int* MmortonCodesIdxUnsortl; //0 1 2 3 ... nbodies-1		
	int* MmortonCodesIdxl;

	int* MlevelUnsortl;
	int* MlevelSortl;

	int* MindexUnsortl;  //0 1 2 3 ... nbodies-2
	int* MindexSortl;
	int* MindexSortTl;


	realPoint* Mposl;  //Положения внутренних узлов в дерева Карраса
	realPoint* Msizel; //Размеры внутренних ячеек
	int* Mparentl;     //Номер ячейки-родителя
	intPair* Mchildl;  //Потомки внутренних ячеек
	intPair* Mrangel;  //Диапазон частиц во внутренней ячейке

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
		vel.resize(nbodies);

		cft.resize(order * (order + 1));
		setBinom(cft);
	}
	catch (...)
	{
		fprintf(stderr, "cannot allocate host memory\n");  exit(-1);
	}

	vtx.resize(nbodies);

	for (int i = 0; i < nbodies; i++) {
		infile >> vtx[i].r()[0] >> vtx[i].r()[1] >> vtx[i].g();
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
			unsigned long long int mem = 0;
						
			massl = (int*)cudaNew(nbodies - 1, sizeof(int));
				mem += (nbodies - 1) * sizeof(int);

			vtxl = (realVortex*)cudaNew(nbodies, sizeof(realVortex));
				mem += nbodies * sizeof(realVortex);

			momsl = (realPoint*)cudaNew((nbodies-1) * order, sizeof(realPoint));
				mem += (nbodies - 1) * order * sizeof(realPoint);

			

			vell = (realPoint*)cudaNew(nbodies, sizeof(realPoint));
				mem += nbodies * sizeof(realPoint);

			vellBS = (realPoint*)cudaNew(nbodies, sizeof(realPoint));
				mem += nbodies * sizeof(realPoint);

			maxrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
			minrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
				mem += 2 * blocks * FACTOR1 * sizeof(realPoint);

			///For MortonTree
			MmortonCodesKeyUnsortl = (int*)cudaNew(nbodies, sizeof(int));
			MmortonCodesKeyl = (int*)cudaNew(nbodies, sizeof(int));
			MmortonCodesIdxUnsortl = (int*)cudaNew(nbodies, sizeof(int));			
			MmortonCodesIdxl = (int*)cudaNew(nbodies, sizeof(int));
				mem += 4 * nbodies * sizeof(int);

			Mposl = (realPoint*)cudaNew(nbodies - 1, sizeof(realPoint));
			Msizel = (realPoint*)cudaNew(nbodies - 1, sizeof(realPoint));
				mem += 2 * (nbodies - 1) * sizeof(realPoint);

			Mparentl = (int*)cudaNew(nnodes, sizeof(int));
				mem += nnodes * sizeof(int);

			Mchildl = (intPair*)cudaNew(nbodies-1, sizeof(intPair));
				mem += (nbodies - 1) * sizeof(intPair);

			Mrangel = (intPair*)cudaNew(nnodes, sizeof(intPair)); //Нужно ли для всех?
				mem += nnodes * sizeof(intPair);

			MlevelUnsortl = (int*)cudaNew(nbodies - 1, sizeof(int));
			MlevelSortl = (int*)cudaNew(nbodies - 1, sizeof(int));
			MindexUnsortl = (int*)cudaNew(nbodies - 1, sizeof(int));
			MindexSortl = (int*)cudaNew(nbodies - 1, sizeof(int));
			MindexSortTl = (int*)cudaNew(nbodies - 1, sizeof(int));
				mem += 5 * (nbodies - 1) * sizeof(int);

				printf("Total allocated memory: %llu bytes = %llu Mbytes\n", mem, mem >> 20);
		}

		if (run == 0)
		{
			cudaCopyVecToDevice(vtx.data(), vtxl, nbodies, sizeof(realVortex));
			//cudaCopyVecToDevice(vel.data(), vell, nbodies, sizeof(realPoint));
			//cudaCopyVecToDevice(vel.data(), vellBS, nbodies, sizeof(realPoint));
			setBinomCftConst(cft.data());
		}

		// run timesteps (launch GPU kernels)

		double starttime, endtime;
		starttime = omp_get_wtime();


		timing[0] += cuInitializationKernel();
	
		/// \brief Построение габаритного прямоугольника
		/// 
		/// \param[in] nbodies количество вихрей
		/// \param[in] vtxl указатель на массив на device, где хранятся вихри
		/// \param[out] Mposl указатель на массив на device, куда записываются координаты центров внутренних ячеек, заполняется только нулевая ячейка(корень)
		/// \param[out] maxrl указатель на пару чисел на device, куда записываются координаты правого верхнего угла габаритного прямоугольника
		/// \param[out] minrl указатель на пару чисел на device, куда записываются координаты левого нижнего угла габаритного прямоугольника
		/// 
		/// \return время исполнения
		timing[1] += McuBoundingBoxKernel(nbodies, vtxl, Mposl, maxrl, minrl);

		//realPoint maxrh, minrh;
		//cudaCopyVecFromDevice(maxrl, &maxrh, 2, sizeof(real));
		//cudaCopyVecFromDevice(minrl, &minrh, 2, sizeof(real));
		//
		//std::cout << "x : " << minrh[0] << "..." << maxrh[0] << "\n";
		//std::cout << "y : " << minrh[1] << "..." << maxrh[1] << "\n";

		/// \brief Вычисление кодов Мортона
		/// 
		/// \param[in] nbodies количество вихрей
		/// \param[in] vtxl указатель на массив на device, где хранятся вихри
		/// \param[out] MmortonCodesKeyUnsortl указатель на массив на device, куда записываются несортированные коды Мортона (в том же порядке, что и вихри в массиве vtxl)
		/// \param[out] MmortonCodesIdxUnsortl указатель на массив на device, куда записываются числа по порядку от 0 до nbodies-1
		/// \param[out] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
		/// \param[out] MmortonCodesIdxl указатель на массив на device, куда записываются правила перестановки кодов Мортона (получается "синхронной" сортировкой MmortonCodesKeyl) 
		/// \param[out] Mrangel указатель на массив на device, куда записываются диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках (заполняется только для корня)
		/// 
		/// \return время исполнения
		auto ttt1 = McuMortonCodesKernel(nbodies, vtxl, MmortonCodesKeyUnsortl, MmortonCodesIdxUnsortl, MmortonCodesKeyl, MmortonCodesIdxl, Mrangel);
		timing[2] += ttt1;
//		std::cout << "McuMortonCodesKernel = " << ttt1 << std::endl;

		/// \brief Определение топологии Мортоновского дерева 
		/// 
		/// \param[in] nbodies количество вихрей
		/// \param[in] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
		/// \param[out] Mparentl указатель на массив на device, куда записываются индексы родителей ячеек 
		/// \param[out] Mchildl указатель на массив на device, куда записываются пары индексов ячеек потомков 
		/// \param[out] Mrangel указатель на массив на device, куда записываются диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках 
		/// 
		/// \return время исполнения
		auto ttt2 = McuMortonInternalNodesKernel(nbodies, MmortonCodesKeyl, Mparentl, Mchildl, Mrangel);
		timing[2] += ttt2;
//		std::cout << "McuMortonInternalNodesKernel = " << ttt2 << std::endl;

		/// \brief Определение геометрических параметров внутренних ячеек Мортоновского дерева 
		/// 
		/// \param[in] nbodies количество вихрей
		/// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
		/// \param[in] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
		/// \param[out] Mposl указатель на массив на device, куда записываются координаты центров внутренних ячеек
		/// \param[out] Msizel указатель на массив на device, куда записываются пары чисел - размеры внутренних ячеек по горизонтали и вертикали
		/// \param[in] Mrangel указатель на массив на device, где хранятся диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках 
		/// \param[out] MlevelUnsortl указатель на массив на device, куда записываются уровни внутренних ячеек мортоновского дерева
		/// \param[out] MlevelSortl указатель на массив на device, куда записываются отсортированные уровни внутренних ячеек мортоновского дерева
		/// \param[out] MindexUnsortl указатель на массив на device, куда записываются числа по порядку от 0 до nbodies-2
		/// \param[out] MindexSortl указатель на массив на device, куда записываются правила перестановки MlevelSortl (получается "синхронной" сортировкой MlevelSortl) 
		/// \param[out] MindexSortTl указатель на массив на device, куда записываются правила обратной перестановки MlevelSortl
		/// 
		/// \return время исполнения
		auto ttt3 = McuMortonInternalCellsGeometryKernel(nbodies, nnodes, MmortonCodesKeyl, Mposl, Msizel, Mrangel, MlevelUnsortl, MlevelSortl, MindexUnsortl, MindexSortl, MindexSortTl);
		timing[2] += ttt3;
//		std::cout << "McuMortonInternalCellsGeometryKernel = " << ttt3 << std::endl;
		
		/// \brief Обнуление необходимых параметров 
		/// 
		/// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
		/// \param[in] nbodies количество вихрей
		/// \param[out] massl указатель на массив на device, где присваиваются -1 внутренним узлам и 1 вихрям  
		/// \param[out] momsl указатель на массив на device, где обнуляются мультипольные моменты всех ячеек
		/// 
		/// \return время исполнения
		auto ttt4 = cuClearKernel2(nnodes, nbodies, massl, momsl);
		timing[3] += ttt4;
//		std::cout << "McuCleearKernel23 = " << ttt4 << std::endl;

		/// \brief Вычисление мультипольных моментов 
		/// 
		/// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
		/// \param[in] nbodies количество вихрей
		/// \param[in] Mchildl указатель на массив на device, куда записываются пары индексов ячеек потомков 
		/// \param[out] massl указатель на массив на device, где присваиваются массы (единица для точечного вихря, число вихрей для ячейки)
		/// \param[out] momsl указатель на массив на device, куда записываются мультипольные моменты всех ячеек
		/// \param[in] vtxl указатель на массив на device, где хранятся вихри
		/// \param[in] MmortonCodesIdxl указатель на массив на device, где хранятся правила перестановки кодов Мортона 
		/// \param[in] Mposl указатель на массив на device, где хранятся координаты центров внутренних ячеек
		/// \param[in] MindexSortl указатель на массив на device, где хранятся правила перестановки массива уровней дерева (MlevelSortl)
		/// \param[in] MindexSortTl указатель на массив на device, где хранятся правила обратной перестановки массива уровней дерева (MlevelSortl)
		/// 
		/// \return время исполнения
		timing[4] += cuSummarizationKernel2(nnodes, nbodies, Mchildl, massl, momsl, vtxl, MmortonCodesIdxl, Mposl, MindexSortl, MindexSortTl);
		
		/// \brief Вычисление скоростей
		/// 
		/// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
		/// \param[in] nbodies количество вихрей
		/// \param[in] itolsq параметр близости
		/// \param[in] epssq радиус вихря
		/// \param[in] Mchildl указатель на массив на device, где хранятся пары индексов ячеек потомков 
		/// \param[in] momsl указатель на массив на device, где хранятся мультипольные моменты всех ячеек
		/// \param[in] vtxl указатель на массив на device, где хранятся вихри
		/// \param[in] MmortonCodesIdxl указатель на массив на device, где хранятся правила перестановки кодов Мортона 
		/// \param[in] Mposl указатель на массив на device, где хранятся координаты центров внутренних ячеек
		/// \param[in] MindexSortl указатель на массив на device, где хранятся правила перестановки массива уровней дерева (MlevelSortl)
		/// \param[in] MindexSortTl указатель на массив на device, где хранятся правила обратной перестановки массива уровней дерева (MlevelSortl)
		/// \param[out] vell указатель на массив на device, куда записываются скорости вихрей
		/// \param[in] Msizel указатель на массив на device, где хранятся пары чисел - размеры внутренних ячеек по горизонтали и вертикали
		/// 
		/// \return время исполнения
		timing[5] += cuForceCalculationKernel2(nnodes, nbodies, itolsq, epssq, Mchildl, momsl, vtxl, MmortonCodesIdxl, Mposl, MindexSortl, MindexSortTl, vell, Msizel);
				       
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

		PrintStatistics(run, runs, timing, mintiming, avtiming, runtime, minruntime, avruntime);
    }

	if (compare || save)
	{
		cudaCopyVecFromDevice(vell, vel.data(), nbodies, sizeof(realPoint));
		for (int i = 0; i < nbodies; ++i)
			vel[i] *= IDPI;
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
			float tBS = cuForceDirectCalculationKernel(nnodes, nbodies, epssq, vtxl, vellBS); 
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
			absVel = std::max(absVel, veloBS[i].length<real>());
			//err += (vel[i] - veloBS[i]).length2();
		}
		err /= nbodies;

		//err = sqrt(err)/nbodies;
		//absVel = 1.0;
		
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
    
 	cudaDelete(massl);
	cudaDelete(vtxl);
	cudaDelete(vell);
	cudaDelete(vellBS);
    cudaDelete(momsl);
   
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