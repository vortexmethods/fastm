/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Tree.cpp                                                         |
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
\brief Основные операции работы с деревом
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/

#include <algorithm>
#include <tuple>

#include <omp.h>

#include "Tree.h"

namespace BH
{
	extern long long op;

	//Конструктор
	MortonTree::MortonTree(std::vector<PointsCopy>& points)
		: pointsCopy(points)
	{
		mortonTree.resize(2 * points.size());
		mortonCodes.resize(points.size());
		mortonLowCells.reserve(points.size());
		
		mortonCodes_temp.resize(points.size());
		s.resize(256 * omp_get_max_threads());


		//Вычисляем факториалы и биномиальные коэффиценты
		iFact[0] = 1.0;
		for (int i = 1; i <= order; ++i)
			iFact[i] = iFact[i - 1] / i;
		getBinom(0, 0) = 1.0;
		for (int i = 1; i <= order; ++i) {
			getBinom(i, 0) = 1.0;
			getBinom(i, i) = 1.0;
			for (int j = 1; j < i; ++j)
				getBinom(i, j) = getBinom(i - 1, j - 1) + getBinom(i - 1, j);
		}//for i
	}//MortonTree(...)



	//Поиск габаритного прямоугольника системы точек
	std::pair<Point2D, Point2D> MortonTree::find_enclosing_rectangle_old(const std::vector<PointsCopy>& points)
	{
		double shared_minx = 1e+10, shared_maxx = -1e+10, shared_miny = 1e+10, shared_maxy = -1e+10;
#pragma omp parallel 
		{
			double minx = 1e+10, maxx = -1e+10, miny = 1e+10, maxy = -1e+10;
#pragma omp for nowait
			for (int ii = 0; ii < (int)points.size(); ++ii)
			{
				const Point2D& r = points[ii].r();
				minx = std::min(r[0], minx);
				maxx = std::max(r[0], maxx);
				miny = std::min(r[1], miny);				
				maxy = std::max(r[1], maxy);
			}//for ii
#pragma omp critical 
			{
				shared_minx = std::min(shared_minx, minx);
				shared_maxx = std::max(shared_maxx, maxx);
				shared_miny = std::min(shared_miny, miny);
				shared_maxy = std::max(shared_maxy, maxy);
			}//critical
		}//parallel
		
		return { {shared_minx, shared_miny}, {shared_maxx, shared_maxy} };
	}//find_enclosing_rectangle_old(...)


/*
//Поиск габаритного прямоугольника системы точек
#pragma omp declare reduction(min : struct PointsCopy : \
        omp_out.r()[0] = omp_in.r()[0] > omp_out.r()[0]  ? omp_out.r()[0] : omp_in.r()[0], \
        omp_out.r()[1] = omp_in.r()[1] > omp_out.r()[1]  ? omp_out.r()[1] : omp_in.r()[1] ) \
        initializer( omp_priv = Vortex2D{ { 1e+10, 1e+10 } , 0.0} )

#pragma omp declare reduction(max : struct PointsCopy : \
        omp_out.r()[0] = omp_in.r()[0] < omp_out.r()[0]  ? omp_out.r()[0] : omp_in.r()[0],  \
        omp_out.r()[1] = omp_in.r()[1] < omp_out.r()[1]  ? omp_out.r()[1] : omp_in.r()[1] ) \
        initializer( omp_priv = Vortex2D{ { -1e+10, -1e+10 } , 0.0} )

	
	std::pair<Point2D, Point2D> MortonTree::find_enclosing_rectangle(const std::vector<PointsCopy>& points)
	{
		PointsCopy minp(Vortex2D{ { 1e+10, 1e+10 } , 0.0}), maxp(Vortex2D{ { -1e+10, -1e+10 } , 0.0 });
				
#pragma omp parallel for reduction(min:minp) reduction(max:maxp)
		for (int i = 0; i < (int)points.size(); ++i) {
			if (points[i].r()[0] < minp.r()[0]) minp.r()[0] = points[i].r()[0];
			if (points[i].r()[1] < minp.r()[1]) minp.r()[1] = points[i].r()[1];
			if (points[i].r()[0] > maxp.r()[0]) maxp.r()[0] = points[i].r()[0];
			if (points[i].r()[1] > maxp.r()[1]) maxp.r()[1] = points[i].r()[1];
		}//for i
		return { minp.r(), maxp.r() };		
	}//find_enclosing_rectangle(...)
*/


	//Построение корня дерева и задание его общих параметров
	void MortonTree::makeRootMortonTree()
	{
		auto gabs = find_enclosing_rectangle_old(pointsCopy);

		double iLeft = gabs.first[0];
		double iBottom = gabs.first[1];
		
		double iRight = gabs.second[0];
		double iTop = gabs.second[1];

		Point2D posCentre = { 0.5 * (iLeft + iRight), 0.5 * (iBottom + iTop) };

		double quadSide = std::max(iRight - iLeft, iTop - iBottom);
		iQuadSideVar = 1.0 / (quadSide * (1.0 + 1.0 / ((1 << codeLength) - 1)));
		lowLeft = posCentre - 0.5*quadSide*Point2D{ 1.0, 1.0 };
				
		//iQuadSideVar = 1.0;
		//lowLeft = Point2D{ 0.0, 0.0 };
				
#pragma omp parallel for schedule(dynamic, 1000)
		for (int i = 0; i < (int)pointsCopy.size(); ++i)
		{
			Point2D locR = (pointsCopy[i].r() - lowLeft) * iQuadSideVar;
			unsigned int mortonCode = morton2D(locR);
				
			auto& mc = mortonCodes[i];
			mc.key = mortonCode;
			mc.originNumber = i;
			mc.r = pointsCopy[i].r();
			mc.g = pointsCopy[i].g();
		}//for i

		//RSort_Node3(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size());				
		RSort_Parallel(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size(), s.data());

		mortonTree[0].range = { 0, (int)pointsCopy.size() - 1 };
		mortonTree[0].particle = false;
	}//makeRootMortonTree()


	//Функция вычисления длины общей части префиксов двух(а значит - и диапазона) частиц
	int MortonTree::delta(int i, int j)
	{
		if ((j < 0) || (j > pointsCopy.size() - 1))
			return -1;

		if (i > j)
			std::swap(i, j);

		//if ((i < 0) || (j > n-1))
		//    exit(111);

		const unsigned int& ki = mortonCodes[i].key;
		const unsigned int& kj = mortonCodes[j].key;
		
		//Поиск номера самого старшего ненулевого бита в числе c 
		int count = 0;
		for (unsigned int c = (ki ^ kj); c; c >>= 1, ++count);

		if ((!count) && (i != j))
		{
			int addCount = 0;
			//единички к номерам i и j добавлены для совместимости с Wolfram Mathematica, 
			//для кода они не важны, но дерево без них почти наверняка построится по-другому        
			for (unsigned int add = ((i + 1) ^ (j + 1)); add; add >>= 1, ++addCount);
			return 2 * codeLength + (2 * codeLength - addCount);
		}//if ((!count) && (i != j))

		return (2 * codeLength - count);
	}//delta(...)


	//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
	int MortonTree::prefixLength(int cell)
	{
		const auto& range = mortonTree[cell].range;
		return std::min(delta(range[0], range[1]), 2 * codeLength);
	}//prefixLength(...)

	
	//Функция вычисления общего префикса двух частиц
	std::pair<unsigned int, int> MortonTree::prefix(int cell)
	{
		int length = prefixLength(cell);
		unsigned int el = mortonCodes[mortonTree[cell].range[0]].key;
		return { el >> (2 * codeLength - length), length };
	}//prefix(...)


	//Функция вычисления геометрических параметров внутренней ячейки дерева
	void MortonTree::setCellGeometry(int cell)
	{
		int prLength;
		unsigned int pr;
		std::tie<unsigned int, int>(pr, prLength) = prefix(cell);

		Point2D sz = { 1.0 / (double)(1 << ceilhalf(prLength)), 1.0 / (1 << (prLength / 2)) };

		Point2D pos = 0.5 * sz;

		for (int i = 0; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[0] += 1.0 / (1 << (1 + i / 2));

		for (int i = 1; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[1] += 1.0 / (1 << (1 + i / 2));

		mortonTree[cell].centre = pos * (1.0 / iQuadSideVar) + lowLeft;
		mortonTree[cell].size = sz * (1.0 / iQuadSideVar);
		mortonTree[cell].level = prLength;

		mortonTree[cell].leaf = (prLength >= NumOfLevels);
	}//setCellGeometry(...)


	//Функция построения i-й внутренней ячейки дерева
	void MortonTree::buildInternalTreeCell(int i)
	{
		const int n = (int)pointsCopy.size();
		int d = sign(delta(i, i + 1) - delta(i, i - 1));
		int delta_min = delta(i, i - d);
		
		int Lmax = 2;
		while (delta(i, i + Lmax * d) > delta_min)
			Lmax *= 2;

		int L = 0;
		for (int t = (Lmax >> 1); t >= 1; t >>= 1)		
			if (delta(i, i + (L + t) * d) > delta_min)
				L += t;
		
		int j = i + L * d;
		int delta_node = delta(i, j);
		int s = 0;
		for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
		{
			int dl = delta(i, i + (s + t) * d);
			if (dl > delta_node)
				s += t;
		}//for p
		int gamma = i + s * d + std::min(d, 0);

		auto min_max = std::minmax(i, j);

		const int& left = gamma;
		const int& right = gamma + 1;

		treeCellT& treeCell = mortonTree[i];

		// Левый потомок - лист или внутренний узел
		bool ifLeftParticle = (min_max.first == gamma);
		treeCell.child[0] = ifLeftParticle ? n + left : left;
		mortonTree[treeCell.child[0]].range = { min_max.first, gamma };
		mortonTree[treeCell.child[0]].particle = ifLeftParticle;
		mortonTree[treeCell.child[0]].parent = i;

		// Правый потомок - лист или внутренний узел
		bool ifRightParticle = (min_max.second == gamma + 1);
		treeCell.child[1] = ifRightParticle ? n + right : right;
		mortonTree[treeCell.child[1]].range = { gamma + 1, min_max.second };
		mortonTree[treeCell.child[1]].particle = ifRightParticle;
		mortonTree[treeCell.child[1]].parent = i;
	}//buildInternalTreeCell(...)


	//Построение внутренних ячеек дерева
	void MortonTree::buildMortonInternalTree()
	{
		//auto t1 = -omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < pointsCopy.size() - 1; ++q)
			buildInternalTreeCell(q);
		//t1 += omp_get_wtime();

		mortonTree[pointsCopy.size() - 1].leaf = false;
		mortonTree[pointsCopy.size() - 1].particle = false;

		//auto t2 = -omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < pointsCopy.size() - 1; ++q)
			setCellGeometry(q);
		//t2 += omp_get_wtime();

		//std::cout << t1 << " " << t2 << std::endl;
	}//buildMortonInternalTree()


	//Построение верхушек дерева --- отдельных частиц
	void MortonTree::buildMortonParticlesTree()
	{
#pragma omp parallel for		
		for (int q = 0; q < pointsCopy.size(); ++q)
		{
			//auto& pt = pointsCopy[mortonCodes[q].originNumber];
			auto& pt = mortonCodes[q];

			auto& cell = mortonTree[pointsCopy.size() + q];
			cell.centre = pt.r;

			cell.leaf = true;
			cell.size = { 0.0, 0.0 };
		}//for q
	}//buildMortonParticlesTree(...)


	//Заполнение списка нижних вершин: 
	//рекурсивный алгоритм, быстро работает в последовательном варианте		
	void MortonTree::fillMortonLowCells(int cell)
	{
		//std::cout << "cell = " << cell << (mortonTree[cell].leaf ? " leaf" : "") << std::endl;
		if (mortonTree[cell].leaf)
			mortonLowCells.push_back(cell);
		else
		{
			fillMortonLowCells(mortonTree[cell].child[0]);
			fillMortonLowCells(mortonTree[cell].child[1]);
		}//else
	}//fillMortonLowCells(...)


	//Заполнение списка нижних вершин: 
	//нерекурсивный алгоритм, хорошо распараллеливается, но неэффективен в последовательном варианте
	void MortonTree::fillMortonLowCellsA()
	{
#pragma omp parallel
		{
			std::vector<int> good_matches_private;
#pragma omp for nowait
			for (int cell = 0; cell < mortonTree.size(); ++cell)			
				if (!mortonTree[mortonTree[cell].parent].leaf && mortonTree[cell].leaf)				
					good_matches_private.push_back(cell);				
#pragma omp critical
			mortonLowCells.insert(mortonLowCells.end(), good_matches_private.begin(), good_matches_private.end());
		}//parallel
	}//fillMortonLowCellsA()


	/// Вычисление параметров дерева (циркуляций и высших моментов)
	void MortonTree::calculateMortonTreeParams(int cell, int level)
	{
		//if (cell < 0 || cell >= mortonTree.size())
		//	std::cout << "cell error, cell = " << cell << std::endl;

		auto& cl = mortonTree[cell];
		//std::cout << "cell = " << cell << (mortonTree[cell].particle ? " particle" : "") << std::endl;
		if (!cl.particle)
		{
#pragma omp parallel default(none) shared(cl, level) num_threads(2) if (level < maxLevelOmp + 1)
			{				 
#pragma omp for 
				for (int s = 0; s < 2; ++s)
					calculateMortonTreeParams(cl.child[s], level+1);

#pragma omp single
				{
					cl.mom.toZero({0.0, 0.0});

					numvector<Point2D, order> h;
					for (int s = 0; s < 2; ++s)
					{
						auto& chld = mortonTree[cl.child[s]];
						const Point2D& g = chld.mom[0];
						cl.mom[0] += g;						
						h[0] = (chld.centre) - (cl.centre);
						for (int i = 1; i < order; ++i)
						{
							h[i] = multz(h[i - 1], h[0]);
#ifdef calcOp
							op += 4;
#endif
						}//for i

						for (int p = 1; p <= order; ++p) {
							Point2D shiftMom = chld.mom[p];
							for (int k = 1; k <= p; ++k)
							{
								shiftMom += getBinom(p, k) * multz(chld.mom[p - k], h[k - 1]);
#ifdef calcOp
								op += 6;
#endif
							}//for k
							cl.mom[p] += shiftMom;
						}//for m
					}//for s
				}//single
			}//parallel
		}//if(!particle)
		else
		{
#ifdef pan
			CalcPanelsParams();
#else
			//CalcPointsParams();
			//C учетом того, что дерево строится до частиц --- только монопольный момент отличен от нуля
			cl.mom[0][0] = mortonCodes[cl.range[0]].g; // { pointsCopy[mortonCodes[cl.range[0]].originNumber].g(), 0.0 };
			cl.mom[0][1] = 0.0;
			for (int p = 1; p <= order; ++p)
				cl.mom[p].toZero();			
#endif
		}//else
	}//CalculateTreeParams()


	//Основная функция вычисления вихревого влияния
	void MortonTree::influenceComputation(std::vector<Point2D>& result)
	{
#ifndef calcOp
#pragma omp parallel for schedule(dynamic, 100)
#endif
		for (int i = 0; i < (int)mortonLowCells.size(); ++i)
		{
			auto& lowCell = mortonTree[mortonLowCells[i]];
			lowCell.E.toZero({ 0.0, 0.0 });;

			//tree->CalcLocalCoeffToLowLevel(lowTree);
			calcLocalCoeffToLowLevel(mortonLowCells[i], 0, true);

			//switch (type)
			//{
			//case 0:
				//lowTree->CalcConvVeloByBiotSavart();
				calcVeloBiotSavart(mortonLowCells[i]);
			//	break;
			//case 1:
			//	lowTree->CalcInfluenceFromPanels();
			//	break;
			//case 2:
			//	lowTree->CalcInfluenceFromPanelsToPoints();
			//	break;
			//default:
			//	break;
			//}

			//lowTree->CalcInfluenceFromVortexFarCells();
			//calcVeloTaylorExpansion(mortonLowCells[i]);
		}//for i

		int n = (int)pointsCopy.size();

#pragma omp parallel for 
		for (int i = 0; i < n; ++i)
		{
			result[i] = pointsCopy[i].veloCopy;
#ifdef linScheme
			result[i + n] = pointsCopy[i].veloCopyLin;
#endif
		}//for i
	}//influenceComputation(...)



	//Расчет коэффициентов разложения в ряд Тейлора внутри ячейки нижнего уровня
	void MortonTree::calcLocalCoeffToLowLevel(int lowCell, int fromWho, bool calcCloseTrees)
	{
		
		auto& lt = mortonTree[lowCell];

		if (lowCell != fromWho)
		{
			auto& wh = mortonTree[fromWho];

			double h, h0;
			h = wh.size[0] + wh.size[1];
			h0 = lt.size[0] + lt.size[1];

			Point2D& r0 = lt.centre;
			Point2D& r1 = wh.centre;

			Point2D rr = r0 - r1;

			double crit2 = rr.length2();
#ifdef calcOp
			op += 3;
			op += 2;
#endif
			// если выполнен критерий дальности => считаем коэффициенты
			if ((crit2 >= sqr((h0 + h + 2.0 * eps) / theta)))
			{
				Point2D  rr = r0 - r1;
				rr *= 1.0 / rr.length2();

				//double inrm2 = 1.0 / rr.length2();
#ifdef calcOp
				op += 3;
#endif			
				
				Point2D varTheta = /*inrm2 **/ rr;
				for (int diag = 0; diag < order; ++diag)
				{
					for (int q = 0; q <= diag; ++q)
						lt.E[q] += ((q % 2 ? -1.0 : 1.0) * iFact[diag - q]) * multzA(varTheta, wh.mom[diag - q]);
					varTheta = ((diag + 1) /** inrm2*/) * multz(varTheta, rr);
				}
				lt.E[0] += iFact[order] * multzA(varTheta, wh.mom[order]);

			}//if crit2
			else // если не выполнен критерий, то рекурсия 
			{
				if (!wh.leaf)
				{
					calcLocalCoeffToLowLevel(lowCell, wh.child[0], calcCloseTrees);
					calcLocalCoeffToLowLevel(lowCell, wh.child[1], calcCloseTrees);
				}
				else if (calcCloseTrees)
					lt.closeCells.push_back(fromWho);
			}//else if
		}//if (lowTree != this)
		else if (calcCloseTrees)
			lt.closeCells.push_back(lowCell); //себя тоже добавляем в ближнюю зону 
	}//CalcConvCoeffToLowLevel(...)



	//Расчет влияния от ближних ячеек по формуле Био - Савара
	void MortonTree::calcVeloBiotSavart(int lowCell)
	{
		auto& lc = mortonTree[lowCell];
		for (size_t k = 0; k < lc.closeCells.size(); ++k)
		{
			//Локальные переменные для цикла
			Point2D velI;
			double dst2eps;

			for (size_t i = lc.range[0]; i <= lc.range[1]; ++i)
			{
				//PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];

				velI.toZero();

				const Point2D& posI = mortonCodes[i].r; //itDop.r();

				auto& rg = mortonTree[lc.closeCells[k]].range;
				for (size_t j = rg[0]; j <= rg[1]; ++j)
				{
					auto& pt = mortonCodes[j];//pointsCopy[mortonCodes[j].originNumber];
					const Point2D& posJ = pt.r;
					const double& gamJ = pt.g;
					dst2eps = std::max((posI - posJ).length2(), eps2);
					velI += (gamJ / dst2eps) * (posI - posJ).kcross();
#ifdef calcOp
					op += 3;
#endif
				}//for j

				pointsCopy[mortonCodes[i].originNumber].veloCopy += velI;
				//itDop.veloCopy += velI;
			}//for i 
		}//for k
	}//calcVeloBiotSavart(...)


	// Расчет конвективных скоростей внутри дерева нижнего уровня
	void MortonTree::calcVeloTaylorExpansion(int lowCell)
	{
		auto& lc = mortonTree[lowCell];
		for (size_t i = lc.range[0]; i <= lc.range[1]; ++i)
		{
			//PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];
			//Point2D deltaPos = itDop.r() - lc.centre;
			Point2D deltaPos = mortonCodes[i].r - lc.centre;
			Point2D v = lc.E[0];

			Point2D hh = deltaPos;
			for (int q = 1; q < order - 1; ++q) {
				v += iFact[q] * multzA(lc.E[q], hh);
				hh = multz(hh, deltaPos);
#ifdef calcOp
				op += 10;
#endif
			}
			v += iFact[order-1] * multzA(lc.E[order-1], hh);
#ifdef calcOp
			op += 6;
#endif

			pointsCopy[mortonCodes[i].originNumber].veloCopy += Point2D({ -v[1], v[0] });
			pointsCopy[mortonCodes[i].originNumber].veloCopy *= IDPI;
			
#ifdef linScheme
			itDop.veloCopyLin += Point2D({ -vL[1], vL[0] });
			itDop.veloCopyLin *= IDPI;
#endif
		}//for i
	}//CalcConvVeloToLowLevel()
	


	//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
	unsigned int MortonTree::expandBits(unsigned int v)
	{
		// вставит 1 нуль
		v = (v | (v << 8)) & 0x00FF00FF;      //  00000000`00000000`abcdefgh`ijklmnop 
		//                                      | 00000000`abcdefgh`ijklmnop`00000000
		//                                      = 00000000`abcdefgh`XXXXXXXX`ijklmnop
		//                                      & 00000000`11111111`00000000`11111111
		//                                      = 00000000`abcdefgh`00000000`ijklmnop

		v = (v | (v << 4)) & 0x0F0F0F0F;      //  00000000`abcdefgh`00000000`ijklmnop 
		//                                      | 0000abcd`efgh0000`0000ijkl`mnop0000
		//                                      = 0000abcd`XXXXefgh`0000ijkl`XXXXmnop
		//                                      & 00001111`00001111`00001111`00001111
		//                                      = 0000abcd`0000efgh`0000ijkl`0000mnop

		v = (v | (v << 2)) & 0x33333333;      //  0000abcd`0000efgh`0000ijkl`0000mnop 
		//                                      | 00abcd00`00efgh00`00ijkl00`00mnop00
		//                                      = 00abXXcd`00efXXgh`00ijXXkl`00mnXXop
		//                                      & 00110011`00110011`00110011`00110011
		//                                      = 00ab00cd`00ef00gh`00ij00kl`00mn00op

		v = (v | (v << 1)) & 0x55555555;      //  00ab00cd`00ef00gh`00ij00kl`00mn00op 
		//                                      | 0ab00cd0`0ef00gh0`0ij00kl0`0mn00op0
		//                                      = 0aXb0cXd`0eXf0gXh`0iXj0kXl`0mXn0oXp
		//                                      & 01010101`01010101`01010101`01010101
		//                                      = 0a0b0c0d`0e0f0g0h`0i0j0k0l`0m0n0o0p
		return v;
	}


	//Мортоновский код для пары из чисел типа double
	//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
	unsigned int MortonTree::morton2D(const Point2D& r)
	{
		/*
		if (x < 0) std::cout << "x*... < 0" << std::endl;
		if (y < 0) std::cout << "y*... < 0" << std::endl;

		if (x * twoPowCodeLength > twoPowCodeLengthMinus)
			std::cout << "x*... > ..." << std::endl;

		if (y * twoPowCodeLength > twoPowCodeLengthMinus)
			std::cout << "y*... > ..." << std::endl;

		x = std::min(std::max(x * twoPowCodeLength, 0.0), twoPowCodeLengthMinus);
		y = std::min(std::max(y * twoPowCodeLength, 0.0), twoPowCodeLengthMinus);

		unsigned int xx = expandBits((unsigned int)x);
		unsigned int yy = expandBits((unsigned int)y);
		*/
		const Point2D& rscale = twoPowCodeLength * r;
		const unsigned int& xx = expandBits((unsigned int)(rscale[0]));
		const unsigned int& yy = expandBits((unsigned int)(rscale[1]));
		return yy | (xx << 1);
	}



	//========================================================
	void MortonTree::RSort_step3(TParticleCode* source, TParticleCode* dest, unsigned int n, unsigned int* offset, unsigned char sortable_bit)
	{
		//unsigned char* b = (unsigned char*)&source[n].key + sortable_bit;

		for (unsigned int i = 0; i < n; ++i)
		{
			TParticleCode* src = &source[i];
			unsigned int off = (src->key >> (sortable_bit * 8)) & 0xFF;
			dest[offset[off]++] = *src;
		}
	}
	//========================================================
	void MortonTree::RSort_Node3(TParticleCode* m, TParticleCode* m_temp, unsigned int n)
	{
		// Заводим массив корзин
		unsigned int s[sizeof(m->key) * 256] = { 0 };
		// Заполняем массив корзин для всех разрядов
		for (unsigned int i = 0; i < n; ++i)
		{
			unsigned int key = m[i].key;
			for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
				++s[((key >> (digit << 3)) & 0xFF) + (digit << 8)];
		}

		// Пересчитываем смещения для корзин
		for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
		{
			unsigned int off = 0;
			for (unsigned int i = 0; i < 256; i++)
			{
				auto& rs = s[i + (digit << 8)];
				unsigned int value = rs;
				rs = off;
				off += value;
			}
		}

		// Вызов сортировки по битам от младших к старшим (LSD)
		for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
		{
			RSort_step3(m, m_temp, n, &s[digit << 8], digit);
			std::swap(m, m_temp);
			//TParticleCode* temp = m;
			//m = m_temp;
			//m_temp = temp;
		}

		// Если ключ структуры однобайтовый, копируем отсортированное в исходный массив
		if (sizeof(m->key) == 1)
		{
			std::swap(m, m_temp);
			//TParticleCode* temp = m;
			//m = m_temp;
			//m_temp = temp;
			memcpy(m, m_temp, n * sizeof(TParticleCode));
		}
	}


	void MortonTree::RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s)
	{
		// Количество задействованных потоков
		unsigned char threads = omp_get_max_threads();
		//std::cout << "threads = " << (int)threads << std::endl;

#pragma omp parallel num_threads(threads)
		{
			TParticleCode* source = m;
			TParticleCode* dest = m_temp;
			unsigned int l = omp_get_thread_num();
			unsigned int div = n / omp_get_num_threads();
			unsigned int mod = n % omp_get_num_threads();
			unsigned int left_index = l < mod ? (div + (mod == 0 ? 0 : 1)) * l : n - (omp_get_num_threads() - l) * div;
			unsigned int right_index = left_index + div - (mod > l ? 0 : 1);

			for (unsigned int digit = 0; digit < sizeof(m->key); ++digit)
			{
				unsigned int s_sum[256] = { 0 };
				unsigned int s0[256] = { 0 };
				unsigned char* b1 = (unsigned char*)&source[right_index].key;
				unsigned char* b2 = (unsigned char*)&source[left_index].key;
				while (b1 >= b2)
				{
					++s0[*(b1 + digit)];
					b1 -= sizeof(TParticleCode);
				}
				for (unsigned int i = 0; i < 256; i++)
				{
					s[i + 256 * l] = s0[i];
				}

#pragma omp barrier
				for (unsigned int j = 0; j < threads; j++)
				{
					for (unsigned int i = 0; i < 256; i++)
					{
						s_sum[i] += s[i + 256 * j];
						if (j < l)
						{
							s0[i] += s[i + 256 * j];
						}
					}
				}

				for (unsigned int i = 1; i < 256; ++i)
				{
					s_sum[i] += s_sum[i - 1];
					s0[i] += s_sum[i - 1];
				}
				unsigned char* b = (unsigned char*)&source[right_index].key + digit;
				TParticleCode* v1 = &source[right_index];
				TParticleCode* v2 = &source[left_index];
				while (v1 >= v2)
				{
					dest[--s0[*b]] = *v1--;
					b -= sizeof(TParticleCode);
				}
#pragma omp barrier
				std::swap(source, dest);
			}
		}

		// Если ключ структуры однобайтовый, просто копируем в исходный массив
		if (sizeof(m->key) == 1)
		{
			memcpy(m, m_temp, n * sizeof(TParticleCode));
		}
	}

}//namespace BH

