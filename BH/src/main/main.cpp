/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.2    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/10/22     |
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
\version 1.2
\date 22 октября 2022 г.
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

#include <functional>
#include <variant>

using fA = std::function<void(std::vector<Point2D>&, const std::vector<Vortex2D>&)>;
using fB = std::function<void(std::vector<Point2D>&, const std::vector<Vortex2D>&, const std::vector<Vortex2D>&, const std::vector <Point2D>&, const std::vector<double>&, const std::vector<Vortex2D>& wakeVP)>;

using varfunc = std::variant<fA, fB>;


void CalcVortexVelo();
void CalcVPVelo();
void SolveLinearSystem();

//для расчета wake на wake
void BiotSavartWakeToWake(std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake);
//для расчета wake на wakeVP
void BiotSavartWakeToWakeVP(std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake, const std::vector<Vortex2D>& wakeAirfoil, const std::vector <Point2D>& dataAirfoil, const std::vector<double>& sec, const std::vector<Vortex2D>& wakeVP);

long long BH::op;


void CheckDefines()
{
	int counter = 0;
#ifdef CALCVORTEXVELO
	++counter;
#endif

#ifdef CALCSHEET
	++counter;
#endif

#ifdef CALCVP
	++counter;
#endif	

	if (counter != 1)
	{
		printf("Check defines for type of the problem!\n");
		exit(-1);
	}
}

int main(int argc, char** argv)
{
	CheckDefines();

	PrintLogoToStream(std::cout);
	PrintProperties();
	
#ifdef CALCVORTEXVELO
	//Расчет скоростей, генерируемых облаком точек, в этих же точках
	CalcVortexVelo();
#endif

#ifdef CALCSHEET
	//Решение системы линейных уравнений
	SolveLinearSystem();
#endif

#ifdef CALCVP
	CalcVPVelo();
#endif	
}

void ReadWake(const std::string fileWake, std::vector<Vortex2D>& wakeVortex, bool readGamma)
{
	int n;
	ifstream infile(fileWake);

	infile >> n;
	wakeVortex.resize(n);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < 2; ++j)
			infile >> wakeVortex[i].r()[j];
		if (readGamma)
			infile >> wakeVortex[i].g();
	} //for i
	infile.close();
}//ReadWake(...)


void ReadPanels(const std::string fileAirfoil, std::vector<Vortex2D>& wakeAirfoil, std::vector<Point2D>& dataAirfoil)
{
	//загрузка файла с параметрами панелей
	ifstream infileAirfoil(fileAirfoil);

	int nAirfoil;
	infileAirfoil >> nAirfoil;
	//PrintConfiguration(n);

	wakeAirfoil.resize(nAirfoil);
	dataAirfoil.resize(nAirfoil);
	std::vector<Point2D> tau(nAirfoil);
	std::vector <double> len(nAirfoil);

	for (int i = 0; i < nAirfoil; ++i)
		for (int j = 0; j < 2; ++j)
			infileAirfoil >> dataAirfoil[i][j];

	for (int i = 0; i < nAirfoil; ++i)
	{
		wakeAirfoil[i].r() = 0.5 * (dataAirfoil[(i + 1) % nAirfoil] + dataAirfoil[i]);
		wakeAirfoil[i].g() = 0.0;
	}
	infileAirfoil.close();

	for (int i = 0; i < nAirfoil; ++i)
	{
		tau[i] = dataAirfoil[(i + 1) % nAirfoil] - dataAirfoil[i];
		len[i] = tau[i].normalize();
	}
	infileAirfoil.close();
}//ReadPanels(...)


void ReadPanelsWithSolution(const std::string fileAirfoil, std::vector<Vortex2D>& wakeAirfoil, std::vector<Point2D>& dataAirfoil, std::vector<double>& sec)
{
	ReadPanels(fileAirfoil, wakeAirfoil, dataAirfoil);

#ifndef asympScheme
	ifstream infileSolution(fileAirfoil + "Sol");
#else
	ifstream infileSolution(fileAirfoil + "SolAs");
#endif

	int nAirfoil = wakeAirfoil.size();
	for (int i = 0; i < nAirfoil; ++i)
		infileSolution >> wakeAirfoil[i].g();

#ifdef linScheme
	sec.resize(nAirfoil);
	for (int i = 0; i < nAirfoil; ++i)
		infileSolution >> sec[i];
#endif 

	infileSolution.close();
}//ReadPanelsWithSolution(...)


struct CallBiotSavart {
	const std::vector<Vortex2D>& wake;
	const std::vector<Vortex2D>& wakeAirfoil;
	const std::vector<Point2D>& dataAirfoil;
	const std::vector<double>& sec;
	const std::vector<Vortex2D>& wakeVP;
		
	//std::vector<Point2D>& velo;
	std::vector<Point2D>& veloBS;
	   
	void operator()(fA foo) { foo(veloBS, wake); }
	void operator()(fB foo) { foo(veloBS, wake, wakeAirfoil, dataAirfoil, sec, wakeVP); };

	static const std::vector<Vortex2D>& emptyVortex2D;
	static const std::vector<Point2D>& emptyPoint2D;
	static const std::vector<double>& emptyDouble;

//#pragma warning( push )
#pragma warning( disable : 3295 )
	CallBiotSavart(std::vector<Point2D>& veloBS_, const std::vector<Vortex2D>& wake_) 
		: wake(wake_), veloBS(veloBS_), wakeAirfoil(emptyVortex2D), dataAirfoil(emptyPoint2D), sec(emptyDouble), wakeVP(emptyVortex2D) {};
//#pragma warning( pop )

	CallBiotSavart(std::vector<Point2D>& veloBS_, const std::vector<Vortex2D>& wake_, const std::vector<Vortex2D>& wakeAirfoil_, const std::vector<Point2D>& dataAirfoil_, const std::vector<double>& sec_, const std::vector<Vortex2D>& wakeVP_) 
		: wake(wake_), veloBS(veloBS_), wakeAirfoil(wakeAirfoil_), dataAirfoil(dataAirfoil_), sec(sec_), wakeVP(wakeVP_) {};
};

const std::vector<Vortex2D>& CallBiotSavart::emptyVortex2D = {};
const std::vector<Point2D>& CallBiotSavart::emptyPoint2D = {};
const std::vector<double>& CallBiotSavart::emptyDouble = {};


void CompareVelocity(varfunc BiotSavartFunction, CallBiotSavart& callFunc, const std::vector<Point2D>& velo)
{
	PrintAccuracyHead();

	callFunc.veloBS.resize(velo.size());

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
		std::visit(callFunc, BiotSavartFunction);
		double tBSFinish = omp_get_wtime();

		std::cout << "done!" << std::endl;
		std::cout << "Time (Biot-Savart): " << tBSFinish - tBSStart << " s" << std::endl;

		std::ofstream velFileBS("../../res/velBS" + task + ".txt");
		velFileBS.precision(16);
		for (int i = 0; i < (int)callFunc.veloBS.size(); ++i)
			velFileBS << callFunc.veloBS[i][0] << " " << callFunc.veloBS[i][1] << std::endl;
		velFileBS.close();
	}
	else
	{
		std::ifstream velFileBS("../../res/velBS" + task + ".txt");

		for (int i = 0; i < (int)callFunc.veloBS.size(); ++i)
			velFileBS >> callFunc.veloBS[i][0] >> callFunc.veloBS[i][1];

		velFileBS.close();
		std::cout << "done!" << std::endl;
	}

	double err = 0.0, absVel = 0.0;
	for (int i = 0; i < (int)callFunc.veloBS.size(); ++i)
	{
		err += (callFunc.veloBS[i] - velo[i]).length();
		absVel += callFunc.veloBS[i].length();
	}

	PrintAccuracyError(err / absVel);
}//CompareVelocity


void SaveVelocity(const std::vector<Point2D>& velo)
{
	std::ofstream velFile("../../res/velBH" + task + ".txt");
	velFile.precision(16);
	for (int i = 0; i < (int)velo.size(); ++i)
		velFile << velo[i][0] << " " << velo[i][1] << std::endl;
	velFile.close();
}//SaveVelocity



//Расчет скоростей, генерируемых облаком точек, в этих же точках
void CalcVortexVelo()
{
	std::vector<Vortex2D> wake;
	ReadWake(vortexFile, wake, true); //true - значит, читать тройки (x, y, \Gamma)
	PrintConfiguration(wake.size());

	std::vector<Point2D> velo, veloBS;

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime=1e+10, avruntime = 0.0;

	// Проводим расчет скоростей заданное количество раз
	for (int run = 0; run < runs; ++run) 
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		velo.clear();
		velo.resize(wake.size());

		runtime = -omp_get_wtime();

		// Конструктор для вычисления скоростей частиц
		BarnesHut BH(wake);

		// Построение дерева treeVrt
		BH.BuildNecessaryTrees(timing[1]);
			
		// Расчет влияния вихрей на вихри
		BH.InfluenceComputation(velo, timing[2], timing[3]);

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

	// Сравнение решений по прямому и быстрому методу
	if (compare)
	{
		CallBiotSavart CBS(veloBS, wake);		
		CompareVelocity(BiotSavartWakeToWake, CBS, velo);
	}

	if (save)
		SaveVelocity(velo);

	std::cout << "Goodbye! " << std::endl;
}


void CalcVPVelo()
{
	std::vector<Vortex2D> wakeVortex, wakeAirfoil, wakeVP;
	std::vector<Point2D> velo, veloBS, dataAirfoil;
	
	ReadWake(vortexFile, wakeVortex, true); // Считываем список вихрей

	std::vector<double> sec;
	ReadPanelsWithSolution(airfoilFile, wakeAirfoil, dataAirfoil, sec);

	ReadWake(vpFile, wakeVP, false);  // Считываем список точек, в которых вычисляются скорости

	PrintConfiguration(wakeVP.size());

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime = 1e+10, avruntime = 0.0;

	for (int run = 0; run < runs; ++run)
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		velo.clear();
		velo.resize(wakeVP.size());

		runtime = -omp_get_wtime();

		// Конструктор для вычисления скоростей частиц wakeVР, вызванных влиянием wakeVortex
		BarnesHut BH(wakeVortex, wakeAirfoil, dataAirfoil, sec, wakeVP);

		// Построение деревьев treeVrt, treePan, treeVP
		BH.BuildNecessaryTrees(timing[1]);
				
		// Расчет влияния
		BH.InfluenceComputation(velo, timing[2], timing[3]);
		
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

	// Сравнение решений по прямому и быстрому методу
	if (compare)
	{
		CallBiotSavart CBS(veloBS, wakeVortex, wakeAirfoil, dataAirfoil, sec, wakeVP);
		CompareVelocity(BiotSavartWakeToWakeVP, CBS, velo);				
	}

	if (save)
		SaveVelocity(velo);

	std::cout << "Goodbye! " << std::endl;
}

//Решение системы линейных уравнений
void SolveLinearSystem()
{
	std::vector<Vortex2D> wakeVortex;
	std::vector<Vortex2D> wakeAirfoil;
	std::vector<Point2D> dataAirfoil;
	
	ReadWake(vortexFile, wakeVortex, true);

	ReadPanels(airfoilFile, wakeAirfoil, dataAirfoil);
	int n = wakeAirfoil.size();
	PrintConfiguration(n);

	int vsize = n;
#ifdef  linScheme
	vsize *= 2;
#endif

for (int run = 0; run < runs; ++run)
{
std::vector<double> gam(vsize, 0.0);

	//wake - центры панелей, data - начала и концы панелей
	BarnesHut BH(wakeVortex, wakeAirfoil, dataAirfoil);

	double timing1 = 0.0;

	// Построение деревьев treeVrt, treePan
	BH.BuildNecessaryTrees(timing1);
	
	double timingRhs = 0.0;

	double t1, t2;

	double R = 0.0;	

	std::vector<double> rhs(vsize); 
	std::vector<Point2D> velo;
	velo.resize(vsize);

	// Расчет влияния вихревого следа и набегающего потока на правую часть
	BH.RhsComputation(velo, timing1, timingRhs);
	for (int i = 0; i < n; ++i)
		rhs[i] -= ((velo[i] + velInf) & BH.pointsCopyPan[i].tau);
#ifdef linScheme
	for (int i = 0; i < n; ++i)
		rhs[n+i] -= (velo[n+i] & BH.pointsCopyPan[i].tau);
#endif

	t1 = omp_get_wtime();
	GMRES(BH, gam, R, rhs, n);
	t2 = omp_get_wtime();

	if (run <= 1 )
{
	std::cout << "Tree time = " << timing1 << std::endl;
	std::cout << "GMRES time = " << t2 - t1 << std::endl;
}
	std::cout << "Total time = " << timing1 + t2-t1 << std::endl;

	
	if (compare && run == 0)
	{
#ifdef linScheme
#ifndef asympScheme
		std::ifstream infile("../../res/Exact" + task + "-linear.txt");
#else
		std::ifstream infile("../../res/Exact" + task + "-as.txt");
#endif //asympScheme
#else
		std::ifstream infile("../../res/Exact" + task + "-const.txt");
#endif // linScheme
		

		std::vector<double> gamFile(vsize);
		
		for (size_t i = 0; i < gam.size(); i++) {
			infile >> gamFile[i];
		}
		std::cout << "done!" << std::endl;

		double refVal = 0.0;
		for (int i = 0; i < n; i++) {
		    if (fabs(gam[i]) > refVal)
				refVal = fabs(gam[i]);
		}

#ifndef linScheme
		double err = 0.0, maxErr = 0.0;
		for (int i = 0; i < n; i++) {
			err = fabs(gam[i] - gamFile[i]);
			if (err/refVal > maxErr)
			    maxErr = err/refVal;
		}
		PrintAccuracyError(maxErr);
#else
		double err = 0.0, maxErr = 0.0;
		for (int i = 0; i < n; i++) {
			err = fabs((gam[i] - 0.5*gam[n+i]) - (gamFile[i] - 0.5*gamFile[n+i]));
			if (err/refVal > maxErr)
			    maxErr = err/refVal;
			err = fabs((gam[i] + 0.5*gam[n+i]) - (gamFile[i] + 0.5*gamFile[n+i]));
			if (err/refVal > maxErr)
			    maxErr = err/refVal;
		}
		PrintAccuracyError(maxErr);

#endif
		infile.close();
	}


	if (save && run==0)
	{
		std::ofstream outfile("../../res/gamRes.txt");
		outfile.precision(16);
		for (int i = 0; i < vsize; i++)
			outfile << gam[i] << std::endl;
		outfile << R << std::endl;
		outfile.close();
	}
}
	std::cout << "Goodbye! " << std::endl;
}//SolveLinearSystem()



//для расчета wake на wake
void BiotSavartWakeToWake(std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake)
{
#pragma omp parallel for 
	for (int i = 0; i < (int)wake.size(); ++i)
	{
		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;
		double dst2eps, dst2;
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

		velI *= IDPI;
		velo[i] = velI;
	}//for i 
}//BiotSavartWakeToWake(...)


//для расчета wake на wakeVP
void BiotSavartWakeToWakeVP(std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake, const std::vector<Vortex2D>& wakeAirfoil, const std::vector <Point2D>& dataAirfoil, const std::vector<double>& sec, const std::vector<Vortex2D>& wakeVP)
{
#pragma omp parallel for 
	for (int i = 0; i < (int)wakeVP.size(); ++i)
	{
		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;
		double dst2eps, dst2;

		velI.toZero();

		const Point2D& posI = wakeVP[i].r();

		
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

		
		for (int j = 0; j < (int)wakeAirfoil.size(); ++j)
		{
			const Point2D& posJ = wakeAirfoil[j].r();
			PointsCopy pnl;
			pnl.panBegin = dataAirfoil[j];
			pnl.panEnd = (j < (int)wakeAirfoil.size() - 1) ? dataAirfoil[j + 1] : dataAirfoil[0];
			pnl.len = (pnl.panEnd - pnl.panBegin).length();
			pnl.tau = (1.0 / pnl.len) * (pnl.panEnd - pnl.panBegin);
			pnl.g() = wakeAirfoil[j].g();
#ifdef linScheme			
			pnl.gamLin = sec[j];
#endif
			
			Point2D u0 = pnl.tau;
			Point2D pp = posI - pnl.panEnd;
			Point2D s = posI - pnl.panBegin;
			double alpha = atan2(pp ^ s, pp & s);
			double lambda = 0.5 * log((s & s) / (pp & pp));

			Point2D va = pp + s;
			Point2D vb = pnl.tau;
			Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
			Point2D u1 = (0.5 / pnl.len) * omega;
			
			velI += (pnl.g() / pnl.len) * (-alpha * (u0.kcross()) + lambda * u0).kcross();
#ifdef linScheme						
			velI += (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();

#ifdef asympScheme
			for (int t = 0; t < nAngPoints; ++t)
			{
				// Панель около угловой точки
				if (j == KK[t])
				{
					Point2D u0 = pnl.tau;
					Point2D pp = posI - pnl.panEnd;
					Point2D s = posI - pnl.panBegin;
					double alpha = atan2(pp ^ s, pp & s);
					double lambda = 0.5 * log((s & s) / (pp & pp));

					Point2D va = pp + s;
					Point2D vb = pnl.tau;
					Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
					Point2D u1 = (0.5 / pnl.len) * omega;

					double rotAngleForward = atan2((pnl.panEnd - pnl.panBegin)[1], (pnl.panEnd - pnl.panBegin)[0]);
					Point2D varzFirst = { s[0] * cos(rotAngleForward) + s[1] * sin(rotAngleForward),
						s[1] * cos(rotAngleForward) - s[0] * sin(rotAngleForward) };

					Point2D s1 = sFirst(varzFirst, q[t], p[t], pnl.len);
					Point2D i1as1First = { s1[0] * cos(rotAngleForward) + s1[1] * sin(rotAngleForward),
						-s1[1] * cos(rotAngleForward) + s1[0] * sin(rotAngleForward) };

					Point2D sub = (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();

					Point2D a1 = (-1.0 / (1.0 - (double)p[t] / q[t])) * (-alpha * (u0.kcross()) + lambda * u0).kcross();
					Point2D a2 = i1as1First.kcross();

					Point2D add = (pnl.gamLin / pnl.len) * (a1 + a2);
					velI -= sub;
					velI += add;
				}

				// Панель около угловой точки
				if (j == KKm[t])
				{
					Point2D u0 = pnl.tau;
					Point2D pp = posI - pnl.panEnd;
					Point2D s = posI - pnl.panBegin;
					double alpha = atan2(pp ^ s, pp & s);
					double lambda = 0.5 * log((s & s) / (pp & pp));

					Point2D va = pp + s;
					Point2D vb = pnl.tau;
					Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
					Point2D u1 = (0.5 / pnl.len) * omega;

					double rotAngleBackward = atan2((pnl.panBegin - pnl.panEnd)[1], (pnl.panBegin - pnl.panEnd)[0]);
					Point2D varzLast = { pp[0] * cos(rotAngleBackward) + pp[1] * sin(rotAngleBackward),
						pp[1] * cos(rotAngleBackward) - pp[0] * sin(rotAngleBackward) };

					Point2D s2 = sFirst(varzLast, q[t], p[t], pnl.len);
					Point2D i1as1Last = { s2[0] * cos(rotAngleBackward) + s2[1] * sin(rotAngleBackward),
						-s2[1] * cos(rotAngleBackward) + s2[0] * sin(rotAngleBackward) };

					Point2D sub = (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();
					Point2D add = (pnl.gamLin / pnl.len) * ((-1.0 / (1.0 - (double)p[t] / q[t])) * (-alpha * (u0.kcross()) + lambda * u0).kcross() + i1as1Last.kcross());

					velI -= sub;
					velI += add;
				}
			}
#endif
#endif
		}//for j
		

		velI *= IDPI;
		velo[i] = velI + velInf;
	}//for i 
}//BiotSavart(...)