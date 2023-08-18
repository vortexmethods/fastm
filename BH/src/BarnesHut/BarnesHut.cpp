/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: BarnesHut.cpp                                                    |
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
\brief Основные операции метода Барнса-Хата
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#include <fstream>

#include "BarnesHut.h"

namespace BH
{
	void BarnesHut::CreatePan(const Vortex2D& panCenter, const Point2D& panBegin, const Point2D& panEnd, double gamLin, std::vector<PointsCopy>& pointsCopyPan)
	{
		PointsCopy panel(panCenter);
		
#ifdef linScheme
		panel.gamLin = gamLin;
#endif
		panel.panBegin = panBegin;
		panel.panEnd = panEnd;
		
		panel.tau = panEnd - panBegin;
		panel.len = panel.tau.normalize();
		ADDOP(6);

						
		auto& moms = panel.Mpan;
		moms.resize(prm.order + 1, { 0.0, 0.0 });

		Point2D rcur, rd2Pow, rd2;
		rd2 = 0.5 * (panEnd - panBegin);
		ADDOP(2);

		rcur = rd2Pow = multz(rd2, rd2);
		
		moms[0] = { 1.0, 0.0 };
		for (int k = 2; k <= prm.order; k += 2)
		{
			moms[k] = (1.0 / (k + 1)) * rcur;			
			rcur = multz(rcur, rd2Pow);
			ADDOP(3);
		}		

#ifdef linScheme
		rcur = rd2;
		for (int k = 1; k <= prm.order; k += 2)
		{
			moms[k] = (0.5 / (k + 2)) * rcur;
			rcur = multz(rcur, rd2Pow);
			ADDOP(3);
		}		
#endif

		pointsCopyPan.emplace_back(panel);
	}


	void BarnesHut::ConvertToAsypmPanel(PointsCopy& panel, double mu, bool infinityAtBegin)
	{
#ifdef asympScheme
		auto& moms = panel.MpanAs;
		moms.resize(prm.order + 1, { 0.0, 0.0 });

		double lhalfnum = 1.0;

		for (int num = 1; num <= prm.order; ++num)
		{
			lhalfnum *= (0.5 * panel.len);
			ADDOP(2);

			double sum = 0.0;
			for (int k = 0; k <= num; ++k)
			{
				double prd = 1.0;
				for (int j = 1; j <= k; ++j)
					prd *= 1.0 - (num + 1.0) / j;
				ADDOP(2*k);
				sum += (double)(1 << k) / (mu - 1.0 - k) * prd;
				ADDOP(2);
			}

			if (!(num & 1))
				sum += 1.0 / ((1.0 + num) * (1.0 - mu));
			
			if (!(infinityAtBegin && (num & 1)))
				sum = -sum;

			sum *= lhalfnum;

			double phi = num * atan2(panel.tau[1], panel.tau[0]);
			moms[num] = sum * Point2D{ cos(phi), sin(phi) };
			ADDOP(6);


			if (num & 1)
				panel.Mpan[num].toZero();

			//std::cout << num << " " << moms[num] << std::endl;
		}
#endif
	}


	// Конструктор	для решения интегрального уравнения вместе с правой частью
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<std::vector<Vortex2D>>& pointsPan, const std::vector<std::vector<Point2D>>& panPos)
		: prm(prm_)
	{
		// заполнение массива для панелей
		std::vector<int> n(prm.airfoilFile.size());
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			n[p] = (int)pointsPan[p].size();
		}

		//int n = (int)pointsPan.size();

		pointsCopyPan.resize(prm.airfoilFile.size());
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
			pointsCopyPan[p].reserve(n[p]);
		

		//???

		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			for (int i = 0; i < n[p]; ++i)
				CreatePan(pointsPan[p][i], panPos[p][i], panPos[p][(i + 1) % n[p]], 0.0, pointsCopyPan[p]);

#ifdef linScheme
#ifdef asympScheme
			for (int t = 0; t < prm.nAngPoints; ++t)
			{
				ConvertToAsypmPanel(pointsCopyPan[p][prm.KK[t]], prm.mu[t], true);
				ConvertToAsypmPanel(pointsCopyPan[p][prm.KKm[t]], prm.mu[t], false);
			}
#endif
#endif
		}
		// заполнение массива для вихрей
		pointsCopyVrt.insert(pointsCopyVrt.end(), pointsVrt.begin(), pointsVrt.end());		
	}//BarnesHut(...)



	//Конструктор для вычисления скоростей частиц pointsVrt, влияющих самих на себя
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt)
		: prm(prm_) 
	{
		pointsCopyVrt.insert(pointsCopyVrt.end(), pointsVrt.begin(), pointsVrt.end());
	}//BarnesHut(...)



	//Конструктор для вычисления скоростей частиц pointsVР, вызванных влиянием pointsVrt
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<std::vector<Vortex2D>>& pointsPan, const std::vector<std::vector<Point2D>>& panPos, const std::vector<std::vector<double>>& sec, const std::vector<Vortex2D>& pointsVP)
		: prm(prm_) 
	{
		// заполнение массива для панелей

		std::vector<int> n(prm.airfoilFile.size());
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			n[p] = (int)pointsPan[p].size();
		}

		pointsCopyPan.resize(prm.airfoilFile.size());
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
			pointsCopyPan[p].reserve(n[p]);

		//int n = (int)pointsPan.size();
		//pointsCopyPan.reserve(n);
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			for (int i = 0; i < n[p]; ++i)
			{
#ifdef linScheme
				double linComponent = sec[p][i];
#else
				double linComponent = 0.0;
#endif
				CreatePan(pointsPan[p][i], panPos[p][i], panPos[p][(i + 1) % n[p]], linComponent, pointsCopyPan[p]);
			}

#ifdef linScheme
#ifdef asympScheme
			for (int t = 0; t < prm.nAngPoints; ++t)
			{
				ConvertToAsypmPanel(pointsCopyPan[p][prm.KK[t]], prm.mu[t], true);
				ConvertToAsypmPanel(pointsCopyPan[p][prm.KKm[t]], prm.mu[t], false);
			}
#endif
#endif
		}
		pointsCopyVrt.insert(pointsCopyVrt.end(), pointsVrt.begin(), pointsVrt.end());		
		pointsCopyVP.insert(pointsCopyVP.end(), pointsVP.begin(), pointsVP.end());
	}//BarnesHut(...)


	//Построение одного дерева
	void BarnesHut::BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, bool ifpan, double& time, int index)
	{
		tree = std::make_unique<MortonTree>(prm, maxTreeLevel, pointsCopy, ifpan, index);

		double t1 = omp_get_wtime();
		tree->MakeRootMortonTree();
		tree->BuildMortonInternalTree();
		tree->BuildMortonParticlesTree();

		// Один из двух вариантов заполнения нижних вершин
		//tree->FillMortonLowCellsA(); // Рекурсивный
		tree->FillMortonLowCells();  // Нерекурсивный

		double t2 = omp_get_wtime();

		time += t2-t1;
	}


	//Построение всех нужных деревьев на основе заданных точек pointsCopy  
	void BarnesHut::BuildNecessaryTrees(double& time)
	{
#ifdef needTreeVrt
		BuildOneTree(treeVrt, prm.NumOfLevelsVortex, pointsCopyVrt, false, time, -1); // дерево вихрей
#endif
		
#ifdef needTreePan
		treePan.resize(prm.airfoilFile.size());
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
			BuildOneTree(treePan[p], prm.NumOfLevelsAirfoil, pointsCopyPan[p], true, time, p);  // дерево панелей
#endif
		
#ifdef needTreeVP
		BuildOneTree(treeVP, prm.NumOfLevelsVP, pointsCopyVP, false, time, -2);
#endif
	}
	
	// Обновление циркуляций вихревых элементов (для решения СЛАУ)
	void BarnesHut::UpdateGams(const std::vector<double>& newGam)
	{
		std::vector<int> n(prm.airfoilFile.size());
		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			n[p] = (int)pointsCopyPan[p].size();
		}
		//int n = (int)pointsCopyPan.size();

		int cntr = 0;

		for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
		{
			for (int i = 0; i < n[p]; ++i)
			{
				pointsCopyPan[p][i].g() = newGam[cntr + i];
				pointsCopyPan[p][i].veloCopy.toZero();
#ifdef linScheme
				pointsCopyPan[p][i].gamLin = newGam[cntr + i + n[p]];
				pointsCopyPan[p][i].veloCopyLin.toZero();
#endif
			}
			cntr += n[p];
#ifdef linScheme
			cntr += n[p];
#endif
		}
	}


	// Расчет влияния вихрей на панели
	void BarnesHut::RhsComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl, int p)
	{
			double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
			omp_set_nested(1); 
#else
			// Максимальное число уровней вложенности распараллеливания
			omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif

			auto& treeContr = treePan[p];
			auto& pointsCopy = pointsCopyPan[p];

			treeVrt->CalculateMortonTreeParams(0, 0);
			double tTreeParamsFinish = omp_get_wtime();
			timeParams += tTreeParamsFinish - tTreeParamsStart;
			double tInflStart = omp_get_wtime();


#pragma omp parallel for schedule(dynamic, 10)
			for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
			{
				auto& lci = treeContr->mortonLowCells[i];
				auto& lowCell = treeContr->mortonTree[lci];

				for (auto& e : lowCell.E)
					e.toZero();

				treeContr->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);
				treeContr->CalcVeloBiotSavart(lci, treeVrt);
				treeContr->CalcVeloTaylorExpansion(lci);
				treeContr->mortonTree[lci].closeCells.resize(0);
			}

			double tInflStop = omp_get_wtime();
			timeInfl += tInflStop - tInflStart;

			int n = (int)pointsCopy.size();

#pragma omp parallel for 
			for (int i = 0; i < n; ++i)
			{
				result[i] = IDPI * pointsCopy[i].veloCopy;
				ADDOP(2);
#ifdef CALCSHEET
#ifdef linScheme
				result[i + n] = IDPI * pointsCopy[i].veloCopyLin;
				ADDOP(2);
#endif
#endif
			}//for i
	}



// Расчет влияния 
#ifdef CALCSHEET
void BarnesHut::InfluenceComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl, int pContr, int pInf)
{
	double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
	omp_set_nested(1);
#else
	omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif

	//for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	//{

		auto& treeContr = treePan[pContr];
		auto& pointsCopy = pointsCopyPan[pContr];

		if (pContr == 0)
			treePan[pInf]->CalculateMortonTreeParams(0, 0);

		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();

		for (int ii = 0; ii < (int)pointsCopy.size(); ++ii)
		{
			pointsCopy[ii].veloCopy.toZero();
			pointsCopy[ii].veloCopyLin.toZero();
		}


#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
		{
			auto& lci = treeContr->mortonLowCells[i];
			auto& lowCell = treeContr->mortonTree[lci];

			for (auto& e : lowCell.E)
				e.toZero();
			treeContr->mortonTree[lci].closeCells.resize(0);

			if (treePan[pInf]->index >= 0)
			{
				treeContr->mortonTree[lci].closeCellsPfl.resize(prm.airfoilFile.size());
				//for (int jj = 0; jj < prm.airfoilFile.size(); ++jj)
					treeContr->mortonTree[lci].closeCellsPfl[pInf].resize(0);
			}

			treeContr->CalcLocalCoeffToLowLevel(lci, treePan[pInf], 0, true);

			treeContr->CalcInfluenceFromPanels(lci, treePan[pInf]);

			treeContr->CalcVeloTaylorExpansion(lci);
		}

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = IDPI * pointsCopy[i].veloCopy + prm.velInf;
			ADDOP(2);
#ifdef linScheme
			result[i + n] = IDPI * pointsCopy[i].veloCopyLin;
			ADDOP(2);
#endif
		}//for i
	//}//for p

		}//InfluenceComputation(...)	
#endif

	// Расчет влияния 
#ifndef CALCSHEET
void BarnesHut::InfluenceComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl)
{
	double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
	omp_set_nested(1);
#else
	omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif

#ifdef CALCVORTEXVELO
	auto& treeContr = treeVrt;
	auto& pointsCopy = pointsCopyVrt;
#endif

#ifdef CALCVP
	for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	{
		auto& treeContr = treeVP;
		auto& pointsCopy = pointsCopyVP;
#endif

		treeVrt->CalculateMortonTreeParams(0, 0);

#ifdef CALCVP
		treePan[p]->CalculateMortonTreeParams(0, 0);
#endif

		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();


#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
		{
			auto& lci = treeContr->mortonLowCells[i];
			auto& lowCell = treeContr->mortonTree[lci];

			for (auto& e : lowCell.E)
				e.toZero();

			treeContr->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);
			treeContr->CalcVeloBiotSavart(lci, treeVrt);
			treeContr->CalcVeloTaylorExpansion(lci);



#ifdef CALCVP
			for (auto& e : lowCell.E)
				e.toZero();
			treeContr->mortonTree[lci].closeCells.resize(0);
			treeContr->CalcLocalCoeffToLowLevel(lci, treePan[p], 0, true);

			treeContr->CalcVeloBiotSavart(lci, treePan[p]);

			treeContr->CalcVeloTaylorExpansion(lci);
#endif//defined CALCVP
		}

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 			
		for (int i = 0; i < n; ++i)
		{
			result[i] = IDPI * pointsCopy[i].veloCopy + prm.velInf;
			ADDOP(2);
		}//for i
#ifdef CALCVP
}//for p
#endif

		}//InfluenceComputation(...)
#endif


void BarnesHut::IterativeInfluenceComputation(std::vector<Point2D>& result, const std::vector<double>& newGam, double& timeParams, double& timeInfl, int pContr, int pInf)
{
	//for (size_t p = 0; p < prm.airfoilFile.size(); ++p)
	//{
		//UpdateGams(newGam);

	double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
					omp_set_nested(1);
#else
	omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif

	auto& treeContr = treePan[pContr];
	auto& pointsCopy = pointsCopyPan[pContr];
	if (pContr == 0)
		treePan[pInf]->CalculateMortonTreeParams(0, 0);
	double tTreeParamsFinish = omp_get_wtime();
	timeParams += tTreeParamsFinish - tTreeParamsStart;

	double tInflStart = omp_get_wtime();
	for (int ii = 0; ii < (int)pointsCopy.size(); ++ii)
	{
		pointsCopy[ii].veloCopy.toZero();
		pointsCopy[ii].veloCopyLin.toZero();
	}

#pragma omp parallel for schedule(dynamic, 10)
	for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
	{
		auto& lci = treeContr->mortonLowCells[i];
		auto& lowCell = treeContr->mortonTree[lci];

		for (auto& e : lowCell.E)
			e.toZero();
		
		treeContr->CalcLocalCoeffToLowLevel(lci, treePan[pInf], 0, false);
		treeContr->UpdateInfluence(lci, treePan[pInf]);

		treeContr->CalcVeloTaylorExpansion(lci);
	}

	double tInflStop = omp_get_wtime();
	timeInfl += tInflStop - tInflStart;

	int n = (int)pointsCopyPan[pContr].size();
	#pragma omp parallel for 
	for (int i = 0; i < n; ++i)
	{
		result[i] = IDPI * pointsCopy[i].veloCopy;
		ADDOP(2);
#ifdef linScheme
		result[i + n] = IDPI * pointsCopy[i].veloCopyLin;
		ADDOP(2);
#endif
	}
	//}
}

//*/

}//namespace BH