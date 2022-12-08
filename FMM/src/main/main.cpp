/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: main.cpp                                                         |
| Info: Source code of FMM                                                    |
|                                                                             |
| This file is part of FMM.                                                   |
| FMM is free software: you can redistribute it and/or modify it              |
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
| along with FMM.  If not, see <http://www.gnu.org/licenses/>.                |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Multipole method for 2D vortex particles
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

#include <omp.h>

#include "Cell.h"
#include "FastMultipole.h"
#include "Logo.h"
#include "Particle.h"
#include "Point2D.h"
#include "Tree.h"


using namespace FMM;

const double IDPI = 0.15915494309189534;

int Particle::counter = 0;
void BiotSavart(const std::vector<Particle>& wake, std::vector<Point2D>& velo);

extern double mt, m2mt;

int main(int argc, char** argv)
{
	PrintLogoToStream(std::cout);
	PrintProperties();
	
	std::vector<Particle> pp;

	std::fstream f1(nameFile);
	f1 >> n;
	pp.resize(n);
	for (int i = 0; i < n; ++i)
	{
		f1 >> pp[i].r[0]; f1 >> pp[i].r[1];
		f1 >> pp[i].q;

		pp[i].z = std::complex<double>{ pp[i].r[0], pp[i].r[1] };
	}
	f1.close(); 
	//std::cout << "particles are ready!" << std::endl;

	PrintConfiguration(n);

	double runtime, minruntime = 1e+10, avruntime = 0;
	double timing[7], mintiming[7]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[7]{};

	for (int run = 0; run < runs; ++run)
	{
		FastMultipole::tDn = 0.0;
		FastMultipole::tUp = 0.0;
		FastMultipole::tLeaf = 0.0;
		FastMultipole::tL2L = 0.0;
		FastMultipole::tM2L = 0.0;
		
		for (auto& x : pp)
			x.f = 0;

		double t1 = omp_get_wtime();

		Tree qt(pp);

		//auto maxit = std::max_element(qt.node.begin(), qt.node.end(), [](const FMM::Cell& n1, const Cell& n2) {return n1.level < n2.level; });

		//std::cout << "max_level = " << maxit->level << std::endl;

		double t2 = omp_get_wtime();

		//std::cout << "tree time = " << t2 - t1 << std::endl;

		FastMultipole f(nt);
		f.InfluenceComputation(qt);

		double t3 = omp_get_wtime();

		runtime = t3 - t1;
		if (minruntime > runtime)
			minruntime = runtime;
		avruntime += runtime;

		timing[1] = t2 - t1; //tree
		timing[2] = FastMultipole::tUp; //tUp
		timing[3] = FastMultipole::tM2L; //tM2L
		timing[4] = FastMultipole::tL2L; //tL2L
		timing[5] = FastMultipole::tLeaf; //tLeaf  //tDN = tM2L + tL2L + tLeaf

		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];
		for (int i = 1; i <= 6; ++i)
		{
			if (mintiming[i] > timing[i])
				mintiming[i] = timing[i];
			avtiming[i] += timing[i];
		}
		//std::cout << mt << " " << m2mt << std::endl;
		PrintStatistics(run, runs, timing, mintiming, avtiming, runtime, minruntime, avruntime, f.opDiv + f.opMult + 2LL * n);
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
			BiotSavart(pp, veloBS);
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
			Point2D velo = { pp[i].f.imag() * IDPI, pp[i].f.real() * IDPI };
			err += (velo - veloBS[i]).length();
			absVel += veloBS[i].length();
		}
		
		PrintAccuracyError(err / absVel);
	}

	if (save)
	{
		std::ofstream velFile("../../res/FMM" + task + ".txt");
		velFile.precision(16);
		for (int i = 0; i < n; i++)
			velFile << pp[i].f.imag() * IDPI << " " << pp[i].f.real() * IDPI << std::endl;
		velFile.close();
	}

	std::cout << "Goodbye! " << std::endl;


}



void BiotSavart(const std::vector<Particle>& wake, std::vector<Point2D>& velo)
{
#pragma omp parallel for 
	for (int i = 0; i < (int)wake.size(); ++i)
	{
		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;
		double dst2eps, dst2;
		double cft = IDPI;

		velI.toZero();

		const Point2D& posI = wake[i].r;

		for (int j = 0; j < (int)wake.size(); ++j)
		{
			const Point2D& posJ = wake[j].r;
			const double& gamJ = wake[j].q;

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

