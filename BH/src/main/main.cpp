/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: main.cpp                                                         |
| Info: Source code of BH                                                     |
|                                                                             |
| This file is part of BH.                                                    |
| BH is free software: you can redistribute it and/or modify it               |
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
| along with BH.  If not, see <http://www.gnu.org/licenses/>.                 |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Barnes-Hut method for 2D vortex particles
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/

#include <algorithm>
#include <fstream> 
#include <iostream>
#include <vector>

#include "omp.h"

#include "BarnesHut.h"
#include "Logo.h"
#include "Methods.h"
#include "Params.h"
#include "Tree.h"
#include "Vortex2D.h"

using std::ifstream;
using std::ofstream;

using namespace BH;


void CalcVortexVelo();
void SolveLinearSystem();
void BiotSavart(const std::vector<Vortex2D>& wake, std::vector<Point2D>& velo);

long long BH::op;

int main(int argc, char** argv)
{
	PrintLogoToStream(std::cout);
	PrintProperties();

	CalcVortexVelo();
	//SolveLinearSystem();
}

void CalcVortexVelo()
{
	int n;
	std::vector<Vortex2D> wake;
	std::vector<Point2D> velo;

	ifstream infile;
	infile.open(nameFile);

	infile >> n;

	wake.resize(n);
	//velo.resize(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2; ++j)
			infile >> wake[i].r()[j];
		
		infile >> wake[i].g();
	} //for i

	infile.close();

	PrintConfiguration(n);

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime=1e+10, avruntime = 0.0;

	for (int run = 0; run < runs; ++run)
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		velo.clear();
		velo.resize(n);

		runtime = -omp_get_wtime();

		BarnesHut BH(wake);

		////// 
		BH.BuildTree(timing[1]);
		////// 

		//double timeStartBH, timeStopBH;
		BH.InfluenceComputation(velo, 0, timing[2], timing[3]);


		runtime += omp_get_wtime();
		if (minruntime > runtime)
			minruntime = runtime;
		avruntime += runtime;

		timing[4] = timing[1] + timing[2] + timing[3];
		for (int i = 1; i <= 4; ++i)
		{
			if (mintiming[i] > timing[i])
				mintiming[i] = timing[i];
			avtiming[i] += timing[i];
		}
		
		PrintStatistics(run, runs, timing, mintiming, avtiming, runtime, minruntime, avruntime);
	}

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
		std::ofstream velFile("../../res/BH" + task + ".txt");
		velFile.precision(16);
		for (int i = 0; i < n; i++)
			velFile << velo[i][0] << " " << velo[i][1] << std::endl;
		velFile.close();
	}

	std::cout << "Goodbye! " << std::endl;
}


void BiotSavart(const std::vector<Vortex2D>& wake, std::vector<Point2D>& velo)
{
#pragma omp parallel for 
	for (int i = 0; i < (int)wake.size(); ++i)
	{
		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;
		double dst2eps, dst2;
		double cft = 0.5 / PI;

		velI.toZero();

		const Point2D& posI = wake[i].r();

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



void SolveLinearSystem()
{
	int n;

	ifstream infile;
	infile.open(nameFile);

	infile >> n;
	std::cout << "n = " << n << std::endl;

	std::vector<Vortex2D> wake(n);
	std::vector<Point2D> tau(n), data(n);
	std::vector <double> len(n);
	int vsize = n;
#ifdef  linScheme
	vsize *= 2;
#endif
	std::vector <double> gam(vsize, 0.0);


	for (int i = 0; i < n; ++i)
		for (int j = 0; j < 2; ++j)
			infile >> data[i][j];

	for (int i = 0; i < n - 1; ++i)
	{
		wake[i].r() = 0.5 * (data[i + 1] + data[i]);
		wake[i].g() = gam[i];
	}

	wake[n - 1].r() = 0.5 * (data[0] + data[n - 1]);
	wake[n - 1].g() = gam[n - 1];

	infile.close();

	Point2D tt;
	for (int i = 0; i < n; ++i)
	{
		if (i != n - 1)
			tt = data[i + 1] - data[i];
		else tt = data[0] - data[n - 1];
		len[i] = tt.length();
		tau[i] = (1.0 / len[i]) * tt;
	}

	Point2D velInf = { 0.8660254037844386, -0.5 }; 

	std::vector<double> rhs(vsize);
	for (int i = 0; i < n; ++i)
		rhs[i] = -(velInf & tau[i]);
#ifdef linScheme
	for (int i = n; i < 2 * n; ++i)
		rhs[i] = 0.0;
#endif
	BarnesHut BH(wake, data);
	
	double timing1;
	BH.BuildTree(timing1);

	double t1, t2;

	double R = 0.0;
	t1 = omp_get_wtime();
	//BiCGStab(BH, gam, R, rhs, len, tau, n, 0);
	//GMRES(BH, gam, R, rhs, len, tau, n);
	//GMRESfile(BH, gam, R, rhs, len, tau, 2 * n);
	t2 = omp_get_wtime();

	std::cout << "All Time = " << t2 - t1 << std::endl;

	std::ofstream outfile;
#ifdef pan
	outfile.open("..//res//gamsRes.txt");
#else 
	outfile.open("..//res//gamsRes.txt");
#endif 

	outfile.precision(16);

	for (int i = 0; i < vsize; i++)
		outfile << gam[i] << std::endl;
	outfile << R << std::endl;
	outfile.close();
}

