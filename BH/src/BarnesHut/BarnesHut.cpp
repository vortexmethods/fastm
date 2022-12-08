/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
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
\version 1.3
\date 08 декабря 2022 г.
*/

#include <fstream>

#include "BarnesHut.h"

namespace BH
{
	void BarnesHut::CreatePan(const Vortex2D& panCenter, const Point2D& panBegin, const Point2D& panEnd, double gamLin)
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
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos)
		: prm(prm_)
	{
		// заполнение массива для панелей
		int n = (int)pointsPan.size();
		pointsCopyPan.reserve(n);
		
		for (int i = 0; i < n; ++i)		
			CreatePan(pointsPan[i], panPos[i], panPos[(i + 1) % n], 0.0);

#ifdef linScheme
#ifdef asympScheme
		for (int t = 0; t < prm.nAngPoints; ++t)
		{
			ConvertToAsypmPanel(pointsCopyPan[prm.KK[t]], prm.mu[t], true);		
			ConvertToAsypmPanel(pointsCopyPan[prm.KKm[t]], prm.mu[t], false);
		}
#endif
#endif
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
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos, std::vector<double> sec, const std::vector<Vortex2D>& pointsVP)
		: prm(prm_) 
	{
		// заполнение массива для панелей
		int n = (int)pointsPan.size();
		pointsCopyPan.reserve(n);

		for (int i = 0; i < n; ++i)		
		{
#ifdef linScheme
			double linComponent = sec[i];
#else
			double linComponent = 0.0;
#endif
			CreatePan(pointsPan[i], panPos[i], panPos[(i + 1) % n], linComponent);
		}
		
#ifdef linScheme
#ifdef asympScheme
		for (int t = 0; t < prm.nAngPoints; ++t)
		{
			ConvertToAsypmPanel(pointsCopyPan[prm.KK[t]], prm.mu[t], true);
			ConvertToAsypmPanel(pointsCopyPan[prm.KKm[t]], prm.mu[t], false);
		}
#endif
#endif
		pointsCopyVrt.insert(pointsCopyVrt.end(), pointsVrt.begin(), pointsVrt.end());		
		pointsCopyVP.insert(pointsCopyVP.end(), pointsVP.begin(), pointsVP.end());
	}//BarnesHut(...)


	//Построение одного дерева
	void BarnesHut::BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, bool ifpan, double& time)
	{
		tree = std::make_unique<MortonTree>(prm, maxTreeLevel, pointsCopy, ifpan);

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
		BuildOneTree(treeVrt, prm.NumOfLevelsVortex, pointsCopyVrt, false, time); // дерево вихрей
#endif
		
#ifdef needTreePan
		BuildOneTree(treePan, prm.NumOfLevelsAirfoil, pointsCopyPan, true, time);  // дерево панелей
#endif
		
#ifdef needTreeVP
		BuildOneTree(treeVP, prm.NumOfLevelsVP, pointsCopyVP, false, time);
#endif
	}
	
	// Обновление циркуляций вихревых элементов (для решения СЛАУ)
	void BarnesHut::UpdateGams(std::vector<double>& newGam)
	{
		int n = (int)pointsCopyPan.size();

		for (int i = 0; i < n; ++i)
		{
			pointsCopyPan[i].g() = newGam[i];
			pointsCopyPan[i].veloCopy.toZero();
#ifdef linScheme
			pointsCopyPan[i].gamLin = newGam[i + n];
			pointsCopyPan[i].veloCopyLin.toZero();
#endif
		}
	}


	// Расчет влияния вихрей на панели
	void BarnesHut::RhsComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl)
	{
		double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
		omp_set_nested(1);
#else
		// Максимальное число уровней вложенности распараллеливания
		omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif

		auto& treeContr = treePan;
		auto& pointsCopy = pointsCopyPan;
		
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

#ifdef CALCSHEET
		auto& treeContr = treePan;
		auto& pointsCopy = pointsCopyPan;
#endif

#ifdef CALCVP
		auto& treeContr = treeVP;
		auto& pointsCopy = pointsCopyVP;
#endif

#if defined CALCVORTEXVELO || defined CALCVP
		treeVrt->CalculateMortonTreeParams(0, 0);
#endif

#if defined CALCSHEET || defined CALCVP
		treePan->CalculateMortonTreeParams(0, 0);
#endif

		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		//tree->tempBuffer.resize(pointsCopy.size());

		

		//std::cout << "size = " << (int)treeContr->mortonLowCells.size() << std::endl;
		double tInflStart = omp_get_wtime();


#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
		{
			auto& lci = treeContr->mortonLowCells[i];
			auto& lowCell = treeContr->mortonTree[lci];
			


#if defined CALCVORTEXVELO || defined CALCVP
			for (auto& e : lowCell.E)
				e.toZero();
			
			//treeContr->mortonTree[lci].closeCells.resize(0);
			treeContr->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);
			treeContr->CalcVeloBiotSavart(lci, treeVrt);
			treeContr->CalcVeloTaylorExpansion(lci);
#endif //defined CALCVORTEXVELO || defined CALCVP



#if defined CALCSHEET || defined CALCVP
			for (auto& e : lowCell.E)
				e.toZero();
			treeContr->mortonTree[lci].closeCells.resize(0);
			treeContr->CalcLocalCoeffToLowLevel(lci, treePan, 0, true);

#ifdef CALCVP 
			treeContr->CalcVeloBiotSavart(lci, treePan);
#endif
		
#ifdef CALCSHEET
			treeContr->CalcInfluenceFromPanels(lci);
#endif
			
			treeContr->CalcVeloTaylorExpansion(lci);
#endif//defined CALCSHEET || defined CALCVP
		}

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = IDPI * pointsCopy[i].veloCopy + prm.velInf;
			ADDOP(2);
#ifdef CALCSHEET
#ifdef linScheme
			result[i + n] = IDPI * pointsCopy[i].veloCopyLin;
			ADDOP(2);
#endif
#endif
		}//for i
	}//InfluenceComputation(...)


	

	//Расчет влияния вихревого следа на правую часть СЛАУ 
	void BarnesHut::FillRhs(std::vector<double>& rhs, double& timeParams, double& timeInfl)
	{
		double tTreeParamsStart = omp_get_wtime();

		auto& pointsCopy = pointsCopyPan;

		treeVrt->CalculateMortonTreeParams(0, 0);
		double tTreeParamsFinish = omp_get_wtime();

		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();

#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < (int)treePan->mortonLowCells.size(); ++i)
		{
			auto& lci = treePan->mortonLowCells[i];

			auto& lowCell = treePan->mortonTree[lci];
			for (auto& e : lowCell.E)
				e.toZero();

			treePan->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);

			treePan->CalcInfluenceFromPointsToPanel(lci, treeVrt);
			lowCell.closeCells.clear();

			treePan->CalcVeloTaylorExpansion(lci);

			treePan->AddVelo(lci);
		}//for i

		double tInflStop = omp_get_wtime();
		timeInfl = tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			rhs[i] -= pointsCopy[i].velTau;
#ifdef linScheme
			rhs[i + n] -= pointsCopy[i].velTauLin;
#endif
		}//for i
	}//FillRhs(...)


	void BarnesHut::IterativeInfluenceComputation(std::vector<Point2D>& result, std::vector<double>& newGam, double& timeParams, double& timeInfl)
	{
		UpdateGams(newGam);
		
		double tTreeParamsStart = omp_get_wtime();
		
#ifdef OLD_OMP
		omp_set_nested(1);
#else
		omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif
		
		treePan->CalculateMortonTreeParams(0, 0);
		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < (int)treePan->mortonLowCells.size(); ++i)
		{
			auto& lci = treePan->mortonLowCells[i];

			auto& lowCell = treePan->mortonTree[lci];
			for (auto& e : lowCell.E)
				e.toZero();

			treePan->CalcLocalCoeffToLowLevel(lci, treePan, 0, false);
			treePan->UpdateInfluence(lci);

			treePan->CalcVeloTaylorExpansion(lci);
		}

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopyPan.size();
#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = IDPI * pointsCopyPan[i].veloCopy;
			ADDOP(2);
#ifdef linScheme
			result[i + n] = IDPI * pointsCopyPan[i].veloCopyLin;
			ADDOP(2);
#endif
		}
	}
}//namespace BH