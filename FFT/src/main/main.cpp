/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: main.cpp                                                         |
| Info: Source code of FFT                                                    |
|                                                                             |
| This file is part of FFT.                                                   |
| FFT is free software: you can redistribute it and/or modify it              |
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
| along with FFT.  If not, see <http://www.gnu.org/licenses/>.                |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief FFT-based method for 2D vortex particles
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
*/


#include <fstream>
#include <iostream>
#include <vector>

#include "omp.h"

#include "FastFourier.h"
#include "Logo.h"
#include "Params.h"

using namespace FFT;

void BiotSavart(const std::vector<Vortex2D> &wake, std::vector<Point2D>& velo);

int main(int argc, char* argv[])
{
	PrintLogoToStream(std::cout);
	PrintProperties();
	
	int maxTh = omp_get_max_threads();
	   	 
	double 
		tCalcMeshStart, tCalcMeshFinish, \
		tFillGammaStart, tFillGammaFinish, \
		tCorrStart, tCorrFinish;

	double tPuassStart, tPuassFinish;

	double tTotVelStart, tTotVelFinish;

	double tExclude;

	Point2D lowLeft = { 1.0 / (NP - 1), 1.0 / (NP - 1) };
	Point2D dim = { 1 - 2.0 / (NP - 1), 1 - 2.0 / (NP - 1) };
	numvector<int, 2> nNodes = { (NP - 2), (NP - 2) };

	int n;
	std::vector<Vortex2D> wake;

	std::ifstream file(nameFile);
	file >> n;
	wake.resize(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2; j++)
			file >> wake[i].r()[j];
		file >> wake[i].g();
	}
	file.close();

	PrintConfiguration(n);

	double hx = dim[0] / (nNodes[0] - 1);

	wake.resize(n);

	std::vector<Point2D> velo;

	double timing[6], mintiming[6]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[6]{};
	double runtime, minruntime = 1e+10, avruntime = 0.0;

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	for (int run = 0; run < runs; ++run)
	{
		velo.clear();
		velo.resize(n, { 0.0, 0.0 });

		FastFourier fft(wake, velo, lowLeft, dim, nNodes, interpType);
		fft.ReadMatrC();

		omp_set_num_threads(numTh1 ? numTh1 : maxTh);

		tCalcMeshStart = omp_get_wtime();
		fft.InitializationCells();
		tCalcMeshFinish = omp_get_wtime();

		//std::cout << "CalcMesh " << tCalcMeshFinish - tCalcMeshStart << std::endl;

		tFillGammaStart = omp_get_wtime();
		fft.FillGam();
		tFillGammaFinish = omp_get_wtime();

		//std::cout << "FillGamma " << tFillGammaFinish - tFillGammaStart << std::endl;

		tCorrStart = omp_get_wtime();
		fft.CalcCorrectionVel();
		tCorrFinish = omp_get_wtime();

		//std::cout << "VelCorr " << tCorrFinish - tCorrStart << std::endl;

		tExclude = 0.0; //заполнение матрицы Q(r) и ее FFT

		tPuassStart = omp_get_wtime();
		fft.PuassonSolver(tExclude);
		tPuassFinish = omp_get_wtime();

		//std::cout << "FFT " << tPuassFinish - tPuassStart - tExclude << std::endl;

		omp_set_num_threads(numTh3 ? numTh3 : maxTh);
		//std::cout << "omp_Q3_threads = " << omp_get_max_threads() << std::endl;

		tTotVelStart = omp_get_wtime();
		fft.CalcTotalVelo();
		tTotVelFinish = omp_get_wtime();

		//std::cout << "CalcTotalVelo " << tTotVelFinish - tTotVelStart << std::endl;
	
		runtime = -tCalcMeshStart + tCalcMeshFinish \
			    - tFillGammaStart + tFillGammaFinish \
				- tCorrStart + tCorrFinish \
				- tPuassStart + tPuassFinish \
				- tTotVelStart + tTotVelFinish;

		if (minruntime > runtime)
			minruntime = runtime;
		avruntime += runtime;
		
		timing[1] = -tCalcMeshStart + tCalcMeshFinish \
			-tFillGammaStart + tFillGammaFinish;
		timing[2] = -tCorrStart + tCorrFinish;
		timing[3] = -tPuassStart + tPuassFinish - tExclude;
		timing[4] = -tTotVelStart + tTotVelFinish;
		timing[5] = timing[1] + timing[2] + timing[3] + timing[4];
		
		for (int i = 1; i <= 5; ++i)
		{
			if (mintiming[i] > timing[i])
				mintiming[i] = timing[i];
			avtiming[i] += timing[i];
		}

		PrintStatistics(run, runs, timing, mintiming, avtiming, runtime, minruntime, avruntime);

	}

//	double t1 = 
//		- tCalcMeshStart + tCalcMeshFinish \
//		- tFillGammaStart + tFillGammaFinish \
//		- tCorrStart + tCorrFinish;
//	double t2 = -tPuassStart + tPuassFinish - tExclude;
//	double t3 = -tTotVelStart + tTotVelFinish;
//
//	std::cout << "_________________________________ " << std::endl;
//	std::cout << "Q1  = " << t1 << std::endl;
//	std::cout << "FFT = " << t2 << std::endl;
//	std::cout << "Q3  = " << t3 << std::endl;
//	std::cout << "_________________________________ " << std::endl;
//	std::cout << "TIME ALL " << t1 + t2 + t3 << std::endl;



	if (compare)
	{
		PrintAccuracyHead();

		std::vector<Point2D> veloBS(n);
		
		bool realBSfromFile = BSfromFile;
		if (BSfromFile)
		{
			std::ifstream f("../../res/velBS" + task + ".txt");
			if (realBSfromFile = f.good())
				std::cout << "File with Biot-Savart results is found, loading it... ";
			else
				std::cout << "File with Biot-Savart results is not found, computing it... ";			
		}
		
		if (!realBSfromFile)
		{
			omp_set_num_threads(maxTh);
			double tBSStart = omp_get_wtime();
			BiotSavart(wake, veloBS);
			double tBSFinish = omp_get_wtime();
			std::cout << "done!" << std::endl;
			std::cout << "Time (Biot-Savart): " << tBSFinish - tBSStart << " s" << std::endl;

			std::ofstream velFileBS("../../res/velBS" + task + ".txt");
			velFileBS.precision(16);
			for (int i = 0; i < n; i++)
				velFileBS << veloBS[i][0] << " " << veloBS[i][1] << std::endl;
			velFileBS.close();
		}
		else
		{
			std::ifstream velFileBS("../../res/velBS" + task + ".txt");

			for (int i = 0; i < n; i++) {
				velFileBS >> veloBS[i][0] >> veloBS[i][1];
			}
			velFileBS.close();
			std::cout << "done!" << std::endl;
		}

		double err = 0.0, absVel = 0.0;
		for (int i = 0; i < n; i++) {
			err += (velo[i] - veloBS[i]).length();
			absVel += veloBS[i].length();
		}

		PrintAccuracyError(err / absVel);
	}
	
	if (save)
	{
		std::ofstream velFile("../../res/FFT" + task + ".txt");
		velFile.precision(16);
		for (int i = 0; i < n; i++)
			velFile << velo[i][0] << " " << velo[i][1] << std::endl;
		velFile.close();
	}

	std::cout << "Goodbye! " << std::endl;

}//main



void BiotSavart(const std::vector<Vortex2D> &wake, std::vector<Point2D>& velo)
{
#pragma omp parallel for	
	for (int i = 0; i < (int)wake.size(); ++i)
	{
		Point2D velI;
		Point2D tempVel;
		double dst2eps, dst2;
		double cft = 0.5 / PI;

		velI.toZero();

		const Point2D &posI = wake[i].r();

		for (int j = 0; j < (int)wake.size(); ++j)
		{
			const Point2D& posJ = wake[j].r();
			const double& gamJ = wake[j].g();

			tempVel.toZero();

			dst2 = dist2(posI, posJ);
			dst2eps = std::max(dst2, eps2);

			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
		}//for j

		velI *= cft;
		velo[i] = velI;
	}//for i 
}//BiotSavart(...)

