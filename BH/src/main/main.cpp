/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
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
\version 1.5
\date 19 июня 2024 г.
*/

#include "BarnesHut.h"
#include "Logo.h"
#include "Methods.h"
#include "Preprocessor.h"
#include "StreamParser.h"


//#define USE_EIGEN
//#define USE_MKL


#ifdef USE_EIGEN
#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "unsupported/Eigen/IterativeSolvers"
#endif

#ifdef USE_MKL
#include <mkl.h>
#endif

using namespace BH;

#include <functional>
#include <variant>


using fA = std::function<void(const params& prm, std::vector<Point2D>&, const std::vector<Vortex2D>&)>;
using fB = std::function<void(const params& prm, std::vector<Point2D>&, const std::vector<Vortex2D>&, const std::vector<std::vector<Vortex2D>>&, const std::vector<std::vector<Point2D>>&, const std::vector<std::vector<double>>&, const std::vector<Vortex2D>& wakeVP)>;

using varfunc = std::variant<fA, fB>;

void CalcVortexVelo(const params& prm);
void CalcVPVelo(const params& prm);
void SolveLinearSystem(const params& prm);
void SolveLinearSystemDirect(const params& prm);

//для расчета wake на wake
void BiotSavartWakeToWake(const params& prm, std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake);
//для расчета wake на wakeVP
void BiotSavartWakeToWakeVP(const params& prm, std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake, const std::vector<std::vector<Vortex2D>>& wakeAirfoil, const std::vector<std::vector<Point2D>>& dataAirfoil, const std::vector<std::vector<double>>& sec, const std::vector<Vortex2D>& wakeVP);



long long BH::op, op_fast, op_direct;

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

void Start(std::unique_ptr<VMlib::StreamParser>& parser, params& prm)
{
#ifdef calcOp
	omp_set_num_threads(1);
#ifdef OLD_OMP
	omp_set_nested(0);
#else
	// Максимальное число уровней вложенности распараллеливания
	omp_set_max_active_levels(0);
#endif
#endif


	parser->get("vortexFile", prm.vortexFile);
	parser->get("NumOfLevelsVortex", prm.NumOfLevelsVortex);

	parser->get("wakeShift", prm.wakeShift);

#if defined CALCSHEET || defined CALCVP
	parser->get("airfoilFile", prm.airfoilFile);
	parser->get("NumOfLevelsAirfoil", prm.NumOfLevelsAirfoil);

	parser->get("vpFile", prm.vpFile);
	parser->get("NumOfLevelsVP", prm.NumOfLevelsVP);

#ifdef asympScheme
	parser->get("nAngPoints", prm.nAngPoints);
	parser->get("KK", prm.KK);
	parser->get("KKm", prm.KKm);
	parser->get("p", prm.p);
	parser->get("q", prm.q);

	for (int i = 0; i < prm.nAngPoints; ++i)
	{
		prm.mu.push_back((double)prm.p[i] / (double)prm.q[i]);
	}
#endif //asympScheme
#endif 

	parser->get("eps", prm.eps);
	prm.eps2 = prm.eps * prm.eps;
	parser->get("velInf", prm.velInf);
	parser->get("task", prm.task);
	parser->get("order", prm.order);
	parser->get("theta", prm.theta);
	parser->get("compare", prm.compare);
	parser->get("runs", prm.runs);
#ifdef calcOp
	prm.runs = 1;
#endif

	parser->get("BSfromFile", prm.BSfromFile);
	parser->get("save", prm.save);

	int coresForTree;
	parser->get("coresForTree", coresForTree);
	if (coresForTree == 0)
		coresForTree = omp_get_max_threads();
#ifdef calcOp
	coresForTree = 1;
#endif
	prm.maxLevelOmp = (int)std::ceil(log2(coresForTree)) - 1;

#if defined CALCSHEET
	parser->get("epsGMRES", prm.epsGMRES);
	parser->get("residualShow", prm.residualShow);
	parser->get("NumOfLevelsVP", prm.NumOfLevelsVP);
#endif 

	PrintLogoToStream(std::cout);
	PrintProperties(prm);

#ifdef CALCVORTEXVELO	
	CalcVortexVelo(prm);
#endif

#ifdef CALCSHEET	
	//Решение системы линейных уравнений	
	//SolveLinearSystemDirect(prm);
	SolveLinearSystem(prm);
#endif

#ifdef CALCVP		
	CalcVPVelo(prm);
#endif
}


int main(int argc, char** argv)
{
	CheckDefines();

	const std::string optFileDefault = "options.txt";	
	params prm;
	
	VMlib::LogStream info;
	info.assignStream(&std::cout, "info");		
	std::stringstream optionsFile(VMlib::Preprocessor(std::string("../") + ((argc > 1) ? std::string(argv[1]) : optFileDefault)).resultString);
	std::stringstream defaultFile(VMlib::Preprocessor(std::string("../defaults.txt")).resultString);
	auto parser = std::make_unique<VMlib::StreamParser>(info, "parser", optionsFile, defaultFile);

	Start(parser, prm);
}


void ReadWake(const std::string fileWake, std::vector<Vortex2D>& wakeVortex, bool readGamma, const Point2D& shift)
{
	int n;
	std::ifstream infile(fileWake);

	if (!infile.good())
	{
		std::cout << "File " << fileWake << " is not found. Error!" << std::endl;
		exit(-2);
	}
	infile >> n;
	wakeVortex.resize(n);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < 2; ++j)
			infile >> wakeVortex[i].r()[j];
		wakeVortex[i].r() += shift;

		if (readGamma)		
			infile >> wakeVortex[i].g();				
	} //for i
	infile.close();
}//ReadWake(...)


void ReadPanels(const std::string fileAirfoil, std::vector<Vortex2D>& wakeAirfoil, std::vector<Point2D>& dataAirfoil, const Point2D& shift)
{
	//загрузка файла с параметрами панелей
	std::ifstream infileAirfoil(fileAirfoil);
	if (!infileAirfoil.good())
	{
		std::cout << "File " << fileAirfoil << " is not found. Error!" << std::endl;
		exit(-2);
	}

	int nAirfoil;
	infileAirfoil >> nAirfoil;
	//PrintConfiguration(n);

	wakeAirfoil.resize(nAirfoil);
	dataAirfoil.resize(nAirfoil);
	std::vector<Point2D> tau(nAirfoil);
	std::vector <double> len(nAirfoil);

	for (int i = 0; i < nAirfoil; ++i)
	{
		for (int j = 0; j < 2; ++j)
			infileAirfoil >> dataAirfoil[i][j];
		dataAirfoil[i] += shift;
	}

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


void ReadPanelsWithSolution(const params& prm, std::vector<std::vector<Vortex2D>>& wakeAirfoil, std::vector<std::vector<Point2D>>& dataAirfoil, std::vector<std::vector<double>>& sec, const Point2D& shift)
{
	for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	{
		ReadPanels(prm.airfoilFile[p], wakeAirfoil[p], dataAirfoil[p], shift);

#ifdef linScheme
#ifdef asympScheme
		std::string suffix = "-as";
#else
		std::string suffix = "-linear";
#endif//asympScheme
#else
		std::string suffix = "-const";
#endif//linScheme

		std::ifstream infileSolution("../../res/Exact" + prm.task + suffix + ".txt");
		if (!infileSolution.good())
		{
			std::cout << "File " << "../../res/Exact" << prm.task << suffix << ".txt" << " is not found. Error!" << std::endl;
			exit(-2);
		}

		int nAirfoil = (int)wakeAirfoil[p].size();
		for (int i = 0; i < nAirfoil; ++i)
			infileSolution >> wakeAirfoil[p][i].g();

#ifdef linScheme
		for (int k = 0; k < prm.airfoilFile.size(); ++k)
			sec[k].resize(nAirfoil);

		for (int i = 0; i < nAirfoil; ++i)
			infileSolution >> sec[p][i];
#endif 

		infileSolution.close();
	}
}//ReadPanelsWithSolution(...)


struct CallBiotSavart {
	const params& prm;
	const std::vector<Vortex2D>& wake;
	const std::vector<std::vector<Vortex2D>>& wakeAirfoil;
	const std::vector<std::vector<Point2D>>& dataAirfoil;
	const std::vector<std::vector<double>>& sec;
	const std::vector<Vortex2D>& wakeVP;
		
	//std::vector<Point2D>& velo;
	std::vector<Point2D>& veloBS;
	   
	void operator()(fA foo) { foo(prm, veloBS, wake); }
	void operator()(fB foo) { foo(prm, veloBS, wake, wakeAirfoil, dataAirfoil, sec, wakeVP); };

	static const std::vector<std::vector<Vortex2D>>& emptyVortex2Dvec;
	static const std::vector<Vortex2D>& emptyVortex2D;
	static const std::vector<std::vector<Point2D>>& emptyPoint2Dvec;
	static const std::vector<std::vector<double>>& emptyDoublevec;

//#pragma warning( push )
//#pragma warning( disable : 3295 )
	CallBiotSavart(const params& prm_, std::vector<Point2D>& veloBS_, const std::vector<Vortex2D>& wake_)
		: prm(prm_), wake(wake_), veloBS(veloBS_), wakeAirfoil(emptyVortex2Dvec), dataAirfoil(emptyPoint2Dvec), sec(emptyDoublevec), wakeVP(emptyVortex2D) {};
//#pragma warning( pop )

	CallBiotSavart(const params& prm_, std::vector<Point2D>& veloBS_, const std::vector<Vortex2D>& wake_, const std::vector<std::vector<Vortex2D>>& wakeAirfoil_, const std::vector<std::vector<Point2D>>& dataAirfoil_, const std::vector<std::vector<double>>& sec_, const std::vector<Vortex2D>& wakeVP_)
		: prm(prm_), wake(wake_), veloBS(veloBS_), wakeAirfoil(wakeAirfoil_), dataAirfoil(dataAirfoil_), sec(sec_), wakeVP(wakeVP_) {};
};

const std::vector<Vortex2D>& CallBiotSavart::emptyVortex2D = {};
const std::vector<std::vector<Vortex2D>>& CallBiotSavart::emptyVortex2Dvec = {};
const std::vector<std::vector<Point2D>>& CallBiotSavart::emptyPoint2Dvec = {};
const std::vector<std::vector<double>>& CallBiotSavart::emptyDoublevec = {};


void CompareVelocity(const params& prm, varfunc BiotSavartFunction, CallBiotSavart& callFunc, const std::vector<Point2D>& velo, const std::string& suffix)
{
	PrintAccuracyHead();

	callFunc.veloBS.resize(velo.size());

	bool realBSfromFile = prm.BSfromFile;
	if (prm.BSfromFile)
	{
		std::ifstream f("../../res/velBS" + prm.task + suffix + ".txt");
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

		std::ofstream velFileBS("../../res/velBS" + prm.task + suffix + ".txt");
		velFileBS.precision(16);
		for (int i = 0; i < (int)callFunc.veloBS.size(); ++i)
			velFileBS << callFunc.veloBS[i][0] << " " << callFunc.veloBS[i][1] << std::endl;
		velFileBS.close();
	}
	else
	{
		std::ifstream velFileBS("../../res/velBS" + prm.task + suffix + ".txt");

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


void SaveVelocity(const params& prm, const std::vector<Point2D>& velo, const std::string& suffix)
{
	std::ofstream velFile("../../res/velBH" + prm.task + suffix + ".txt");
	velFile.precision(16);
	for (int i = 0; i < (int)velo.size(); ++i)
		velFile << velo[i][0] << " " << velo[i][1] << std::endl;
	velFile.close();
}//SaveVelocity


#ifdef CALCVORTEXVELO
//Расчет скоростей, генерируемых облаком точек, в этих же точках
void CalcVortexVelo(const params& prm)
{
	std::vector<Vortex2D> wake;
	ReadWake(prm.vortexFile, wake, true, prm.wakeShift); //true - значит, читать тройки (x, y, \Gamma)
	PrintConfiguration(prm, (int)wake.size(), { 0 }, 0);

	std::vector<Point2D> velo, veloBS;

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime=1e+10, avruntime = 0.0;

	// Проводим расчет скоростей заданное количество раз
	for (int run = 0; run < prm.runs; ++run)
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		velo.clear();
		velo.resize(wake.size());

		runtime = -omp_get_wtime();

		// Конструктор для вычисления скоростей частиц
		BarnesHut BH(prm, wake);

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

		if (run == 0)
			PrintTreesInfo(BH.treeVrt, BH.treePan, BH.treeVP);

		PrintStatistics(run, prm.runs, timing, mintiming, avtiming, runtime, minruntime, avruntime, 0);
#ifdef calcOp	
		PrintOps(BH::op);
		BH::op = 0;
#endif
	}//run


	std::string suffix = "";
	if (prm.save)
		SaveVelocity(prm, velo, suffix);

	// Сравнение решений по прямому и быстрому методу
	if (prm.compare)
	{
		CallBiotSavart CBS(prm, veloBS, wake);		
		CompareVelocity(prm, BiotSavartWakeToWake, CBS, velo, suffix);
#ifdef calcOp	
		PrintOps(BH::op);
		BH::op = 0;
#endif
	}

	std::cout << "Goodbye! " << std::endl;
}
#endif

#ifdef CALCVP
void CalcVPVelo(const params& prm)
{

	std::vector<Vortex2D> wakeVortex, wakeVP;
	std::vector<std::vector<Vortex2D>> wakeAirfoil;
	std::vector<Point2D> velo, veloBS;
	std::vector<std::vector<Point2D>> dataAirfoil;
	
	wakeAirfoil.resize(prm.airfoilFile.size());
	dataAirfoil.resize(prm.airfoilFile.size());
	
	ReadWake(prm.vortexFile, wakeVortex, true, prm.wakeShift); // Считываем список вихрей

	std::vector<std::vector<double>> sec;
	sec.resize(prm.airfoilFile.size()); //???

	ReadPanelsWithSolution(prm, wakeAirfoil, dataAirfoil, sec);

	ReadWake(prm.vpFile, wakeVP, false, {0.0, 0.0});  // Считываем список точек, в которых вычисляются скорости
	
	std::vector<int> n;
	for (auto& x : wakeAirfoil)
		n.push_back(x.size());

	PrintConfiguration(prm, (int)wakeVortex.size(), n, (int)wakeVP.size());

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime = 1e+10, avruntime = 0.0;

	for (int run = 0; run < prm.runs; ++run)
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		velo.clear();
		velo.resize(wakeVP.size());

		runtime = -omp_get_wtime();

		// Конструктор для вычисления скоростей частиц wakeVР, вызванных влиянием wakeVortex
		BarnesHut BH(prm, wakeVortex, wakeAirfoil, dataAirfoil, sec, wakeVP);

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

		if (run == 0)
			PrintTreesInfo(BH.treeVrt, BH.treePan, BH.treeVP);

		PrintStatistics(run, prm.runs, timing, mintiming, avtiming, runtime, minruntime, avruntime, 0);
#ifdef calcOp	
		PrintOps(BH::op);
		BH::op = 0;
#endif
	}//runs

#ifdef linScheme
#ifdef asympScheme
	std::string suffix = "-as";
#else
	std::string suffix = "-linear";
#endif//asympScheme
#else
	std::string suffix = "-const";
#endif//linScheme

	if (prm.save)
		SaveVelocity(prm, velo, suffix);

	// Сравнение решений по прямому и быстрому методу
	if (prm.compare)
	{
		CallBiotSavart CBS(prm, veloBS, wakeVortex, wakeAirfoil, dataAirfoil, sec, wakeVP);
		CompareVelocity(prm, BiotSavartWakeToWakeVP, CBS, velo, suffix);		
#ifdef calcOp	
		PrintOps(BH::op);
		BH::op = 0;
#endif
	}

	std::cout << "Goodbye! " << std::endl;
}
#endif


#ifdef CALCSHEET
//Решение системы линейных уравнений
void SolveLinearSystem(const params& prm)
{
	std::vector<Vortex2D> wakeVortex;
	std::vector<std::vector<Vortex2D>> wakeAirfoil(prm.airfoilFile.size());
	std::vector<std::vector<Point2D>> dataAirfoil(prm.airfoilFile.size());

	ReadWake(prm.vortexFile, wakeVortex, true, prm.wakeShift);

	std::vector<int> n(prm.airfoilFile.size());

	Point2D shifts[] = {{0.0, 0.0}, {0.0, 1.1}, {1.1, 1.1}, {1.1, 0.0} , {1.1, -1.1}, {0.0, -1.1},/* {0.0, -2.2}, {1.1, -2.2}*/ {-1.1, -1.1}, {-1.1, 0.0}, {-1.1, 1.1}};
	for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	{		
		ReadPanels(prm.airfoilFile[p], wakeAirfoil[p], dataAirfoil[p], shifts[p]);
		n[p] = (int)wakeAirfoil[p].size();
	}


	PrintConfiguration(prm, (int)wakeVortex.size(), n, 0);

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime = 1e+10, avruntime = 0.0;

	std::vector<int> vsize = n;
	int totalVsize = 0;
	for (auto& x : vsize)
		totalVsize += x;

#ifdef  linScheme
	for (auto& x : vsize)
		x *= 2;
#endif

	std::vector<std::vector<double>> gam(prm.airfoilFile.size());
	for (int k = 0; k < prm.airfoilFile.size(); ++k)
		gam[k].resize(vsize[k], 0.0);

	std::vector<double> R(prm.airfoilFile.size(), 0.0);

	for (int run = 0; run < prm.runs; ++run)
	{
		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		R.clear();
		R.resize(prm.airfoilFile.size(), 0.0);

		runtime = -omp_get_wtime();

		//wakeAirfoil - центры панелей, dataAirfoil - начала и концы панелей
		double tt1 = -omp_get_wtime();
		BarnesHut BH(prm, wakeVortex, wakeAirfoil, dataAirfoil);
		tt1 += omp_get_wtime();
		std::cout << "T_BH = " << tt1 << std::endl;



		// Построение деревьев treeVrt, treePan
		double tt2 = -omp_get_wtime();
		BH.BuildNecessaryTrees(timing[1]);
		tt2 += omp_get_wtime();
		std::cout << "T_Build = " << tt2 << std::endl;

		//for (int k = 0; k < prm.airfoilFile.size(); ++k)
		//{
		//	gam[k].resize(0);
		//	gam[k].resize(vsize[k], 0.0);

		//	std::ifstream oldSol("../../res/gamRes0.txt");

		//	//for (int p = 0; p < vsize[k]; ++p)
		//	//	oldSol >> gam[k][p];

		//	oldSol.close();
		//}

		double tt3 = -omp_get_wtime();
		std::vector<double> flattenGam;

		for (int k = 0; k < prm.airfoilFile.size(); ++k)
			flattenGam.insert(flattenGam.end(), gam[k].begin(), gam[k].end());

		BH.UpdateGams(flattenGam);


		std::vector<std::vector<double>> rhs(prm.airfoilFile.size());
		std::vector<std::vector<Point2D>> velo(prm.airfoilFile.size());

		for (int k = 0; k < prm.airfoilFile.size(); ++k)
		{
			rhs[k].resize(vsize[k]);
			velo[k].resize(vsize[k]);
		}

		tt3 += omp_get_wtime();
		std::cout << "T_Prep = " << tt3 << std::endl;

		double tt4 = -omp_get_wtime();
		// Расчет влияния вихревого следа и набегающего потока на правую часть
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
			BH.RhsComputation(velo[p], timing[2], timing[3], (int)p);

		for (int p = 0; p < prm.airfoilFile.size(); ++p) {
			for (int i = 0; i < n[p]; ++i)
				rhs[p][i] -= ((velo[p][i] + prm.velInf) & BH.pointsCopyPan[p][i].tau);
#ifdef linScheme
			for (int i = 0; i < n[p]; ++i)
				rhs[p][n[p] + i] -= (velo[p][n[p] + i] & BH.pointsCopyPan[p][i].tau);
#endif
		}//for p
		tt4 += omp_get_wtime();
		std::cout << "T_Rhs = " << tt4 << std::endl;

		double tt5 = -omp_get_wtime();
		double tt = 0.0;
		int niter;
		GMRES(BH, gam, R, rhs, n, prm, timing[2], timing[3], niter);
		tt5 += omp_get_wtime();
		std::cout << "T_GMRES = " << tt5 << std::endl;

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

		if (run == 0)
			PrintTreesInfo(BH.treeVrt, BH.treePan, BH.treeVP);

		PrintStatistics(run, prm.runs, timing, mintiming, avtiming, runtime, minruntime, avruntime, niter);
#ifdef calcOp	
		PrintOps(BH::op);
		BH::op = 0;
#endif
	}//runs

	if (prm.save)
	{
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			std::ofstream outfile("../../res/gamRes" + std::to_string(p) + ".txt");
			outfile.precision(16);
			for (int i = 0; i < vsize[p]; i++)
				outfile << gam[p][i] << std::endl;

			outfile << R[p] << std::endl;

			outfile.close();
		}
	}//save

	if (prm.compare)
	{
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			PrintAccuracyHead();

#ifdef linScheme
#ifndef asympScheme
			std::ifstream infile("../../res/Exact" + prm.task + "-linear-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-linear.txt");
#else
			std::ifstream infile("../../res/Exact" + prm.task + "-as-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-as.txt");
#endif //asympScheme
#else
			std::ifstream infile("../../res/Exact" + prm.task + "-const-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-const.txt");
#endif // linScheme

			if (infile.good())
				std::cout << "File with exact solution is found, loading it... ";
			else
			{
				std::cout << "File with exact solution is not found, error!" << std::endl;
				exit(-1);
			}

			double refVal = 0.0;
		
			std::vector<std::vector<double>> gamFile(prm.airfoilFile.size());
#ifndef linScheme
			gamFile[p].resize(vsize[p]);
#else
			gamFile[p].resize(2 * vsize[p]);
#endif

			//std::vector<double> gamFile(totalVsize);

				for (size_t i = 0; i < vsize[p]; ++i) 
					infile >> gamFile[p][i];

			std::cout << "done!" << std::endl;

				for (int i = 0; i < n[p]; i++)
					if (fabs(gam[p][i]) > refVal)
						refVal = fabs(gam[p][i]);



#ifndef linScheme
			double err = 0.0, maxErr = 0.0;
			int cntr = 0;

				cntr = 0;
				for (int i = 0; i < n[p]; ++i) {
					//std::cout << "i = " << i << " " << gam[p][i] << " " << gamFile[p][cntr] << std::endl;
					err = fabs(gam[p][i] - gamFile[p][cntr]);
					if (err / refVal > maxErr)
						maxErr = err / refVal;
					++cntr;
				}
			//}
			PrintAccuracyError(maxErr);
#else
			double err = 0.0, maxErr = 0.0;
			int cntr = 0;

				for (int i = 0; i < n[p]; ++i) {
					err = fabs(gam[p][i] - gamFile[p][cntr]);
					err = fabs((gam[p][i] - 0.5 * gam[p][n[p] + i]) - (gamFile[p][cntr] - 0.5 * gamFile[p][vsize[p] / 2 + cntr]));
					if (err / refVal > maxErr)
						maxErr = err / refVal;
					err = fabs((gam[p][i] + 0.5 * gam[p][n[p] + i]) - (gamFile[p][cntr] + 0.5 * gamFile[p][vsize[p] / 2 + cntr]));
					if (err / refVal > maxErr)
						maxErr = err / refVal;
					++cntr;
				}
			PrintAccuracyError(maxErr);

#endif
			infile.close();
		}
	}//compare


	std::cout << "Goodbye! " << std::endl;
}//SolveLinearSystem()




//Решение системы линейных уравнений
void SolveLinearSystemDirect(const params& prm)
{
	std::vector<Vortex2D> wakeVortex;
	std::vector<std::vector<Vortex2D>> wakeAirfoil(prm.airfoilFile.size());
	std::vector<std::vector<Point2D>> dataAirfoil(prm.airfoilFile.size());

	ReadWake(prm.vortexFile, wakeVortex, true, prm.wakeShift);

	std::vector<int> n(prm.airfoilFile.size());

	Point2D shifts[] = { {0.0, 0.0}, {0.0, 1.1}, {1.1, 1.1}, {1.1, 0.0}, {1.1, -1.1}, {0.0, -1.1}, {-1.1, -1.1}, {-1.1, 0.0}/*, {-1.1, 1.1}*/ };
	
	std::vector<std::vector<double>> Len(prm.airfoilFile.size());
	std::vector<std::vector<Point2D>> Tau(prm.airfoilFile.size());
	

	for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	{
		ReadPanels(prm.airfoilFile[p], wakeAirfoil[p], dataAirfoil[p], shifts[p]);
		n[p] = (int)wakeAirfoil[p].size();

		Len[p].resize(n[p]);
		Tau[p].resize(n[p]);

		for (int i = 0; i < n[p]; ++i)
		{
			Tau[p][i] = Point2D(dataAirfoil[p][(i + 1) % n[p]] - dataAirfoil[p][i]);
			Len[p][i] = Tau[p][i].normalize();
		}
	}

	PrintConfiguration(prm, (int)wakeVortex.size(), n, 0);

	double timing[5], mintiming[5]{ 1e+10, 1e+10, 1e+10, 1e+10, 1e+10 }, avtiming[5]{};
	double runtime, minruntime = 1e+10, avruntime = 0.0;

	std::vector<int> vsize = n;
	int totalVsize = 0;
	for (auto& x : vsize)
		totalVsize += x;

#ifdef  linScheme
	for (auto& x : vsize)
		x *= 2;
#endif

	std::vector<std::vector<double>> gam(prm.airfoilFile.size());
	for (int k = 0; k < prm.airfoilFile.size(); ++k)
		gam[k].resize(vsize[k], 0.0);

	std::vector<double> R(prm.airfoilFile.size(), 0.0);

	for (int run = 0; run < prm.runs; ++run)
	{
		R.clear();
		R.resize(prm.airfoilFile.size(), 0.0);

		for (int i = 0; i < 5; ++i)
			timing[i] = 0.0;

		runtime = -omp_get_wtime();

		size_t nAllVars = prm.airfoilFile.size();
		for (auto& x : vsize)
			nAllVars += x;

		auto after = [](size_t i, size_t j, size_t n) -> bool
		{
			return ((i == j + 1) || (i == 0 && j == n - 1));
		};//isAfter(...)


		auto locAlpha = [](const Point2D& p, const Point2D& s) -> double
		{
			return atan2(cross3(p, s), p & s);			
		};

		auto locLambda = [](const Point2D& p, const Point2D& s) -> double
		{
			return 0.5 * log((s & s) / (p & p));
		};

		auto locOmega = [](const Point2D& a, const Point2D& b, const Point2D& c) -> Point2D
		{
			return (a & b) * c + (Point2D({ -c[1], c[0] })) * cross3(a, b);			
		};
		
		
		std::vector<double> mtrDir(nAllVars * nAllVars);		
		int posI = 0, posJ = 0;
		double tt1 = -omp_get_wtime();
		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
		{
			posJ = 0;
			int npI = n[aflI];
			for (int aflJ = 0; aflJ < prm.airfoilFile.size(); ++aflJ)
			{
				int npJ = n[aflJ];
#pragma omp parallel for 
				for (int i = 0; i < n[aflI]; ++i)
				{
					Point2D tauI = Tau[aflI][i];
					double lenI = Len[aflI][i];
					Point2D nrmI = { tauI[1], -tauI[0] };

					//std::cout << "i = " << i << std::endl;
//#pragma omp parallel for
					for (int j = 0; j < n[aflJ]; ++j)
					{
						bool self = (aflI == aflJ);
						if ((i == j) && self)
						{
							mtrDir[(posI + i) * nAllVars + (posJ + j)] = -0.5;

#ifdef linScheme
							mtrDir[(posI + i) * nAllVars + (posJ + npJ + j)] = 0.0;
							mtrDir[(posI + npI + i) * nAllVars + (posJ + j)] = 0.0;
							mtrDir[(posI + npI + i) * nAllVars + (posJ + npJ + j)] = -1.0 / 24.0;
#endif
						}//if i==j
						else
						{
							Point2D tauJ = Tau[aflJ][j];
							double lenJ = Len[aflJ][j];
							Point2D nrmJ = { tauJ[1], -tauJ[0] };

							Point2D p1 = dataAirfoil[aflI][(i + 1) % n[aflI]] - dataAirfoil[aflJ][(j + 1) % n[aflJ]];
							Point2D s1 = dataAirfoil[aflI][(i + 1) % n[aflI]] - dataAirfoil[aflJ][j];
							Point2D p2 = dataAirfoil[aflI][i] - dataAirfoil[aflJ][(j + 1) % n[aflJ]];
							Point2D s2 = dataAirfoil[aflI][i] - dataAirfoil[aflJ][j];
							Point2D di = dataAirfoil[aflI][(i + 1) % n[aflI]] - dataAirfoil[aflI][i];
							Point2D dj = dataAirfoil[aflJ][(j + 1) % n[aflJ]] - dataAirfoil[aflJ][j];

							numvector<double, 3> alpha = { \
								(self && after(j, i, n[aflI])) ? 0.0 : locAlpha(s2, s1), \
								locAlpha(s2, p1), \
								(self && after(i, j, n[aflI])) ? 0.0 : locAlpha(p1, p2) \
							};

							numvector<double, 3> lambda = { \
								(self && after(j, i, n[aflI])) ? 0.0 : locLambda(s2, s1), \
								locLambda(s2, p1), \
								(self && after(i, j, n[aflI])) ? 0.0 : locLambda(p1, p2) \
							};

							numvector<Point2D, 3> v00, v11;
							numvector<Point2D, 2> v01, v10;

							Point2D i00, i01, i10, i11;

							v00 = {
								locOmega(s1, tauI, tauJ),
								-locOmega(di, tauI, tauJ),
								locOmega(p2, tauI, tauJ)
							};

							i00 = IDPI / lenI * (-(alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2]).kcross() \
								+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]));

							mtrDir[(posI + i) * nAllVars + (posJ + j)] = i00 & nrmI;

#ifdef linScheme
							v01 = {
								0.5 / (dj.length()) * (((p1 + s1) & tauJ) * locOmega(s1, tauI, tauJ) - s1.length2() * tauI),
								-0.5 * di.length() / dj.length() * locOmega(s1 + p2, tauJ, tauJ)
							};

							i01 = IDPI / lenI * (-((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]).kcross() \
								+ ((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * di.length() * tauJ);

							mtrDir[(posI + i) * nAllVars + (posJ + npJ + j)] = i01 & nrmI;


							v10 = {
								-0.5 / di.length() * (((s1 + s2) & tauI) * locOmega(s1, tauI, tauJ) - s1.length2() * tauJ),
								0.5 * dj.length() / di.length() * locOmega(s1 + p2, tauI, tauI)
							};

							i10 = IDPI / lenI * (-((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1]).kcross() \
								+ ((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * dj.length() * tauI);

							mtrDir[(posI + npI + i) * nAllVars + (posJ + j)] = i10 & nrmI;


							v11 = {
								1.0 / (12.0 * di.length() * dj.length()) * (2.0 * (s1 & locOmega(s1 - 3.0 * p2, tauI, tauJ)) * locOmega(s1, tauI, tauJ) - s1.length2() * (s1 - 3.0 * p2)) - 0.25 * locOmega(s1, tauI, tauJ),
								-di.length() / (12.0 * dj.length()) * locOmega(di, tauJ, tauJ),
								-dj.length() / (12.0 * di.length()) * locOmega(dj, tauI, tauI)
							};

							i11 = IDPI / lenI * (-((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]).kcross()\
								+ (lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
								+ 1.0 / 12.0 * (dj.length() * tauI + di.length() * tauJ - 2.0 * locOmega(s1, tauI, tauJ)));

							mtrDir[(posI + npI + i) * nAllVars + (posJ + npJ + j)] = i11 & nrmI;
#endif
						}
					}

					//std::cout << "j finished" << std::endl;

					if (aflI == aflJ)
					{
						mtrDir[(nAllVars - prm.airfoilFile.size() + aflI) * nAllVars + (posI + i)] = lenI;
						mtrDir[(posI + i) * nAllVars + (nAllVars - prm.airfoilFile.size() + aflI)] = 1.0;
					}
				}//for i

				posJ += vsize[aflI];
			}
			posI += vsize[aflI];
		}
		tt1 += omp_get_wtime();


		double tt2 = -omp_get_wtime();
		std::vector<double> rhsDir(nAllVars);
		posI = 0;
		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
		{
			double velI = 0.0;
			double velILin = 0.0;

			for (int i = 0; i < n[aflI]; ++i)
			{
				const Point2D& posI0 = dataAirfoil[aflI][i];
				const Point2D& posI1 = dataAirfoil[aflI][(i + 1) % n[aflI]];
				Point2D tauI = Tau[aflI][i];
				double lenI = Len[aflI][i];
			
				for (size_t it = 0; it != wakeVortex.size(); ++it)
				{
					const Vortex2D& vt = wakeVortex[it];
					const Point2D& posJ = vt.r();
					const double& gamJ = vt.g();

					Point2D s = posJ - posI0;
					Point2D p = posJ - posI1;
					double alpha = locAlpha(p, s);

					velI -= gamJ * alpha;
#ifdef linScheme
					Point2D u1 = 0.5 / lenI * locOmega(p + s, tauI, tauI);					
					double lambda = locLambda(p, s);
					velILin -= gamJ * (alpha * (u1 & tauI) + lambda * (u1 & (-tauI.kcross())));
#endif
				}

				rhsDir[posI + i] = velI - (prm.velInf & tauI);

#ifdef linScheme
				rhsDir[posI + n[aflI] + i] = velILin;
#endif
			}//panels


			posI += vsize[aflI];
		}
		tt2 += omp_get_wtime();

		/*
		std::ofstream matrFile("matrix.txt");
		//matrFile.precision(16);
		for (int p = 0; p < nAllVars; ++p)
		{
			for (int q = 0; q < nAllVars; ++q)
				matrFile << mtrDir[p * nAllVars + q] << " ";
			matrFile << std::endl;
		}
		matrFile.close();

		std::ofstream rhsFile("rhs.txt");
		//matrFile.precision(16);
		for (int p = 0; p < nAllVars; ++p)		
			rhsFile << rhsDir[p] << std::endl;		
		rhsFile.close();
		//*/

		std::cout << "Eigen::nbThreads() = " << Eigen::nbThreads() << std::endl;

		double tt3;
		
#ifdef USE_EIGEN
		Eigen::MatrixXd EigMtr(nAllVars, nAllVars);
		Eigen::VectorXd EigRhs(nAllVars);

		
		for (int p = 0; p < nAllVars; ++p)
		{
			EigRhs(p) = rhsDir[p];
			for (int q = 0; q < nAllVars; ++q)
				EigMtr(p, q) = mtrDir[p * nAllVars + q];
		}

		/*
		tt3 = -omp_get_wtime();
		Eigen::VectorXd EigSol = EigMtr.partialPivLu().solve(EigRhs);		
		tt3 += omp_get_wtime();
		//*/

		//*
		tt3 =- omp_get_wtime();
		Eigen::GMRES<Eigen::MatrixXd> solver(EigMtr);
		Eigen::VectorXd EigSol = solver.solve(EigRhs);
		tt3 += omp_get_wtime();
		
		std::cout << "#iterations:     " << solver.iterations() << std::endl;
		std::cout << "estimated error: " << solver.error() << std::endl;
		//*/

		posI = 0;
		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
		{
			for (int i = 0; i < vsize[aflI]; ++i)
				gam[aflI][i] = EigSol(posI + i) * Len[aflI][i % n[aflI]];
			posI += vsize[aflI];
		}

		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
			R[aflI] = EigSol(EigSol.size() - prm.airfoilFile.size() + aflI);
#endif


#ifdef USE_MKL
		int* ipiv = new int[nAllVars];
		tt3 = -omp_get_wtime();
		int info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nAllVars, nAllVars, mtrDir.data(), nAllVars, ipiv);
		int info2 = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nAllVars, 1, mtrDir.data(), nAllVars, ipiv, rhsDir.data(), 1);
		tt3 += omp_get_wtime();
		delete[] ipiv;

		posI = 0;
		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
		{
			for (int i = 0; i < vsize[aflI]; ++i)
				gam[aflI][i] = rhsDir[posI + i] * Len[aflI][i % n[aflI]];
			posI += vsize[aflI];
		}

		for (int aflI = 0; aflI < prm.airfoilFile.size(); ++aflI)
			R[aflI] = rhsDir[rhsDir.size() - prm.airfoilFile.size() + aflI];
#endif
		



		/*
		std::ofstream solFile("sol.txt");
		//matrFile.precision(16);
		for (int p = 0; p < nAllVars; ++p)
			solFile << EigSol(p) << std::endl;
		solFile.close();
		*/

		timing[1] = tt1;
		timing[2] = tt2;
		timing[3] = tt3;

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

		std::cout << "timeMatr = " << timing[1] << std::endl;
		std::cout << "timeRhs =  " << timing[2] << std::endl; 
		std::cout << "timeSolve = " << timing[3] << std::endl;
		//if (run == 0)
		//	PrintTreesInfo(BH.treeVrt, BH.treePan, BH.treeVP);

		PrintStatistics(run, prm.runs, timing, mintiming, avtiming, runtime, minruntime, avruntime, 0);
		
		//#ifdef calcOp
		//		PrintOps(BH::op);
		//		BH::op = 0;
		//#endif
	}//runs



	if (prm.save)
	{
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			std::ofstream outfile("../../res/gamRes" + std::to_string(p) + ".txt");
			outfile.precision(16);
			for (int i = 0; i < vsize[p]; i++)
				outfile << gam[p][i] << std::endl;

			outfile << R[p] << std::endl;

			outfile.close();
		}
	}//save

	if (prm.compare)
	{
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			PrintAccuracyHead();

#ifdef linScheme
#ifndef asympScheme
			std::ifstream infile("../../res/Exact" + prm.task + "-linear-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-linear.txt");
#else
			std::ifstream infile("../../res/Exact" + prm.task + "-as-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-as.txt");
#endif //asympScheme
#else
			std::ifstream infile("../../res/Exact" + prm.task + "-const-" + std::to_string(p) + ".txt");
			//std::ifstream infile("../../res/Exact" + prm.task + "-const.txt");
#endif // linScheme

			if (infile.good())
				std::cout << "File with exact solution is found, loading it... ";
			else
			{
				std::cout << "File with exact solution is not found, error!" << std::endl;
				exit(-1);
			}

			double refVal = 0.0;

			std::vector<std::vector<double>> gamFile(prm.airfoilFile.size());
#ifndef linScheme
			gamFile[p].resize(vsize[p]);
#else
			gamFile[p].resize(2 * vsize[p]);
#endif

			//std::vector<double> gamFile(totalVsize);

			for (size_t i = 0; i < vsize[p]; ++i)
				infile >> gamFile[p][i];

			std::cout << "done!" << std::endl;

			for (int i = 0; i < n[p]; i++)
				if (fabs(gam[p][i]) > refVal)
					refVal = fabs(gam[p][i]);



#ifndef linScheme
			double err = 0.0, maxErr = 0.0;
			int cntr = 0;

			cntr = 0;
			for (int i = 0; i < n[p]; ++i) {
				//std::cout << "i = " << i << " " << gam[p][i] << " " << gamFile[p][cntr] << std::endl;
				err = fabs(gam[p][i] - gamFile[p][cntr]);
				if (err / refVal > maxErr)
					maxErr = err / refVal;
				++cntr;
			}
			//}
			PrintAccuracyError(maxErr);
#else
			double err = 0.0, maxErr = 0.0;
			int cntr = 0;

			for (int i = 0; i < n[p]; ++i) {
				err = fabs(gam[p][i] - gamFile[p][cntr]);
				err = fabs((gam[p][i] - 0.5 * gam[p][n[p] + i]) - (gamFile[p][cntr] - 0.5 * gamFile[p][vsize[p] / 2 + cntr]));
				if (err / refVal > maxErr)
					maxErr = err / refVal;
				err = fabs((gam[p][i] + 0.5 * gam[p][n[p] + i]) - (gamFile[p][cntr] + 0.5 * gamFile[p][vsize[p] / 2 + cntr]));
				if (err / refVal > maxErr)
					maxErr = err / refVal;
				++cntr;
			}
			PrintAccuracyError(maxErr);

#endif
			infile.close();
		}
	}//compare

	std::cout << "Goodbye! " << std::endl;
}//SolveLinearSystemDirect()


#endif

//для расчета wake на wake
void BiotSavartWakeToWake(const params& prm, std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake)
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
			dst2eps = std::max(dst2, prm.eps2);

			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
			ADDOP(5);
		}//for j

		velI *= IDPI;
		velo[i] = velI;
		ADDOP(2);
	}//for i 
}//BiotSavartWakeToWake(...)


//для расчета wake на wakeVP
void BiotSavartWakeToWakeVP(const params& prm, std::vector<Point2D>& velo, const std::vector<Vortex2D>& wake, const std::vector<std::vector<Vortex2D>>& wakeAirfoil, const std::vector<std::vector<Point2D>>& dataAirfoil, const std::vector<std::vector<double>>& sec, const std::vector<Vortex2D>& wakeVP)
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
			dst2eps = std::max(dst2, prm.eps2);

			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
			ADDOP(5);
		}//for j

		for (size_t p = 0; p < wakeAirfoil.size(); ++p)
		{
			for (int j = 0; j < (int)wakeAirfoil[p].size(); ++j)
			{
				const Point2D& posJ = wakeAirfoil[p][j].r();
				PointsCopy pnl;
				pnl.panBegin = dataAirfoil[p][j];
				pnl.panEnd = (j < (int)wakeAirfoil[p].size() - 1) ? dataAirfoil[p][j + 1] : dataAirfoil[p][0];
				pnl.len = (pnl.panEnd - pnl.panBegin).length();
				pnl.tau = (1.0 / pnl.len) * (pnl.panEnd - pnl.panBegin);
				pnl.g() = wakeAirfoil[p][j].g();
#ifdef linScheme			
				pnl.gamLin = sec[p][j];
#endif
				ADDOP(6);

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
				ADDOP(43);

#ifdef asympScheme
				for (int t = 0; t < prm.nAngPoints; ++t)
				{
					// Панель около угловой точки
					if (j == prm.KK[t])
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

						Point2D s1 = sFirst(varzFirst, prm.q[t], prm.p[t], pnl.len);
						Point2D i1as1First = { s1[0] * cos(rotAngleForward) + s1[1] * sin(rotAngleForward),
							-s1[1] * cos(rotAngleForward) + s1[0] * sin(rotAngleForward) };

						Point2D sub = (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();

						Point2D a1 = (-1.0 / (1.0 - (double)prm.p[t] / prm.q[t])) * (-alpha * (u0.kcross()) + lambda * u0).kcross();
						Point2D a2 = i1as1First.kcross();

						Point2D add = (pnl.gamLin / pnl.len) * (a1 + a2);
						velI -= sub;
						velI += add;
						ADDOP(56);
					}

					// Панель около угловой точки
					if (j == prm.KKm[t])
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

						Point2D s2 = sFirst(varzLast, prm.q[t], prm.p[t], pnl.len);
						Point2D i1as1Last = { s2[0] * cos(rotAngleBackward) + s2[1] * sin(rotAngleBackward),
							-s2[1] * cos(rotAngleBackward) + s2[0] * sin(rotAngleBackward) };

						Point2D sub = (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();
						Point2D add = (pnl.gamLin / pnl.len) * ((-1.0 / (1.0 - (double)prm.p[t] / prm.q[t])) * (-alpha * (u0.kcross()) + lambda * u0).kcross() + i1as1Last.kcross());

						velI -= sub;
						velI += add;
						ADDOP(56);
					}
				}
#endif
#endif
			}//for j
		}//for p

		velI *= IDPI;
		velo[i] = velI + prm.velInf;
		ADDOP(2);
	}//for i 
}//BiotSavart(...)