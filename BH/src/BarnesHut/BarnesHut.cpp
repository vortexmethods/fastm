/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.0    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2021/08/05     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
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
\version 1.0
\date 05 августа 2021 г.
*/

#include <fstream>

#include "omp.h"

#include "BarnesHut.h"

namespace BH
{

	/// Конструктор	
	BarnesHut::BarnesHut(const std::vector<Vortex2D>& points, const std::vector<Point2D>& panPos)
	{
		int n = (int)points.size();
		pointsCopy.reserve(n);
		double len2;

		for (int i = 0; i < n; ++i)
		{
			pointsCopy.emplace_back(points[i]);

			Point2D& tt = pointsCopy[i].tau;

			pointsCopy[i].gamLin = 0.0;
			pointsCopy[i].panBegin = panPos[i];

			if (i != n - 1)
			{
				pointsCopy[i].panEnd = panPos[i + 1];
				tt = panPos[i + 1] - panPos[i];
			}
			else
			{
				pointsCopy[i].panEnd = panPos[0];
				tt = panPos[0] - panPos[i];
			}
			const double& len = pointsCopy[i].len = tt.length();
			tt /= len;
#ifdef linScheme
			pointsCopy[i].dipP = (len / 12.0) * tt;
#ifdef calcOp
			op += 5;
#endif
#endif
			if (order >= 2)
			{
				len2 = len * len;
				Point2D tau2 = multz(tt, tt);
				pointsCopy[i].quaP = (len2 / 24.0) * tau2;
#ifdef calcOp
				op += 4;
#endif
				if (order >= 3)
				{
					Point2D tau3 = multz(tau2, tt);
#ifdef linScheme
					pointsCopy[i].octP = (len * len2 / 320.0) * tau3;
#ifdef calcOp
					op += 4;
#endif
#endif
					if (order >= 4)
					{
						pointsCopy[i].hexP = (len2 * len2 / 640.0) * multz(tau3, tt);
#ifdef calcOp
						op += 4;
#endif
					}
				}
			}
		}
	}//BarnesHut(...)


	BarnesHut::BarnesHut(const std::vector<Vortex2D>& points)
	{
		pointsCopy.reserve(points.size());

		for (int i = 0; i < (int)points.size(); ++i)
		{
			pointsCopy.emplace_back(points[i]);
		}
	}//BarnesHut(...)

	void BarnesHut::BuildTree(double& time)
	{
		tree.reset(new Tree(pointsCopy));

		tree->makeTree(pointsCopy);

		tree->BigTree = tree.get();

		omp_set_nested(1);

		double t1 = omp_get_wtime();

		tree->CreateTree();
		tree->PushbackLowTrees();

		double t2 = omp_get_wtime();
		//std::cout << "CreateTree time " << t2 - t1 << std::endl;
		time += t2 - t1;

	}

	// Обновление циркуляций вихревых элементов (для решения СЛАУ)
	void BarnesHut::UpdateGams(std::vector<double>& newGam)
	{
		int n = (int)pointsCopy.size();
		for (size_t i = 0; i < n; ++i)
		{
			pointsCopy[i].g() = newGam[i];
			pointsCopy[i].veloCopy.toZero();
#ifdef linScheme
			pointsCopy[i].gamLin = newGam[i + n];
			pointsCopy[i].veloCopyLin.toZero();
#endif
		}
	}

	// Обнуление временной статистики
	void BarnesHut::ClearTimestat()
	{
		tCoeffStart = 0.0;
		tCoeffFinish = 0.0;
		tExactStart = 0.0;
		tExactFinish = 0.0;
		tTreeParamsStart = 0.0;
		tTreeParamsFinish = 0.0;
		tOtherStart = 0.0;
		tOtherFinish = 0.0;
	}

	/// Расчет влияния точек result самих на себя
	void BarnesHut::InfluenceComputation(std::vector<Point2D>& result, int type, double& timeParams, double& timeInfl)
	{
		ClearTimestat();
		tTreeParamsStart += omp_get_wtime();
		tree->CalculateTreeParams();
		tTreeParamsFinish += omp_get_wtime();

		timeParams += tTreeParamsFinish - tTreeParamsStart;
		//std::cout << "CalculateTreeParams() = " << tTreeParamsFinish - tTreeParamsStart << std::endl;

		double tInflStart = omp_get_wtime();
#ifndef calcOp
#pragma omp parallel for schedule(dynamic, 10)
#endif
		for (int i = 0; i < (int)tree->lowTrees.size(); ++i)
		{
			auto& lowTree = tree->lowTrees[i];

			lowTree->ClearCoeff();
			tree->CalcLocalCoeffToLowLevel(lowTree);

			switch (type)
			{
			case 0:
				lowTree->CalcConvVeloByBiotSavart();
				break;
			case 1:
				lowTree->CalcInfluenceFromPanels();
				break;
			case 2:
				lowTree->CalcInfluenceFromPanelsToPoints();
				break;
			default:
				break;
			}
			lowTree->CalcInfluenceFromVortexFarCells();
		}
		double tInflFinish = omp_get_wtime();
		timeInfl += tInflFinish - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = pointsCopy[i].veloCopy;
#ifdef linScheme
			result[i + n] = pointsCopy[i].veloCopyLin;
#endif
		}
	}

	void BarnesHut::IterativeInfluenceComputation(std::vector<Point2D>& result, std::vector<double>& newGam, int type)
	{
		tOtherStart += omp_get_wtime();
		UpdateGams(newGam);
		tOtherFinish += omp_get_wtime();

		tTreeParamsStart += omp_get_wtime();
		tree->CalculateTreeParams();
		tTreeParamsFinish += omp_get_wtime();

#pragma omp parallel for 
		for (int i = 0; i < (int)tree->lowTrees.size(); ++i)
		{
			tOtherStart += omp_get_wtime();
			tree->lowTrees[i]->ClearCoeff();
			tOtherFinish += omp_get_wtime();

			tCoeffStart += omp_get_wtime();
			tree->CalcLocalCoeffToLowLevel(tree->lowTrees[i], false);
			tCoeffFinish += omp_get_wtime();

			tExactStart += omp_get_wtime();
			switch (type)
			{
			case 0:
				tree->lowTrees[i]->CalcConvVeloByBiotSavart();
				break;
			case 1:
				tree->lowTrees[i]->UpdateInfluence();
				break;
			case 2:
				tree->lowTrees[i]->CalcInfluenceFromPanelsToPoints();
				break;
			default:
				break;
			}
			tExactFinish += omp_get_wtime();

			tOtherStart += omp_get_wtime();
			tree->lowTrees[i]->CalcInfluenceFromVortexFarCells();
			tOtherFinish += omp_get_wtime();
		}

		int n = (int)pointsCopy.size();
#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = pointsCopy[i].veloCopy;
#ifdef linScheme
			result[i + n] = pointsCopy[i].veloCopyLin;
#endif
		}
	}


}//namespace BH


