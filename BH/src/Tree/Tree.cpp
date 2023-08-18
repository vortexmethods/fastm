/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
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
\version 1.5
\date 19 июня 2024 г.
*/

#include <algorithm>
#include <cstring>

#include "Tree.h"

namespace BH
{
	extern long long op;

	//Конструктор
	MortonTree::MortonTree(const params& prm_, int maxTreeLevel_, std::vector<PointsCopy>& points, bool ifpans, int index_)
		: prm(prm_), pointsCopy(points), pans(ifpans), maxTreeLevel(maxTreeLevel_), index(index_)
	{
		mortonTree.resize(2 * points.size());
		for (auto& c : mortonTree)
		{
			c.mom.resize(prm.order + 1);
			c.E.resize(prm.order);
		}

		mortonCodes.resize(points.size());
		mortonLowCells.reserve(points.size());
		

		//Вычисляем факториалы и биномиальные коэффиценты
		iFact.resize(prm.order + 1);
		binomCft.resize((prm.order + 1) * (prm.order + 1));
		
		iFact[0] = 1.0;
		for (int i = 1; i <= prm.order; ++i)
			iFact[i] = iFact[i - 1] / i;
		getBinom(0, 0) = 1.0;
		for (int i = 1; i <= prm.order; ++i) 
		{
			getBinom(i, 0) = 1.0;
			getBinom(i, i) = 1.0;
			for (int j = 1; j < i; ++j)
				getBinom(i, j) = getBinom(i - 1, j - 1) + getBinom(i - 1, j);
		}//for i
	}//MortonTree(...)



	//Поиск габаритного прямоугольника системы точек
	std::pair<Point2D, Point2D> MortonTree::FindEnclosingRectangleOld(const std::vector<PointsCopy>& points)
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
	}//FindEnclosingRectangleOld(...)


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

	
	std::pair<Point2D, Point2D> MortonTree::FindEnclosingRectangle(const std::vector<PointsCopy>& points)
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
	}//FindEnclosingRectangle(...)
*/


	//Построение корня дерева и задание его общих параметров
	void MortonTree::MakeRootMortonTree()
	{
		auto gabs = FindEnclosingRectangleOld(pointsCopy);

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
			unsigned int mortonCode = Morton2D(locR);
				
			auto& mc = mortonCodes[i];
			mc.key = mortonCode;
			mc.originNumber = i;
			mc.r = pointsCopy[i].r();
			mc.g = pointsCopy[i].g();
		}//for i

		//std::cout << "omp_threads = " << omp_get_max_threads() << std::endl;

		//RSort_Node3(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size());				
		
		//Временные массивы для сортировки		
		std::vector<unsigned int> s(256 * omp_get_max_threads());
		std::vector<TParticleCode> mortonCodes_temp(pointsCopy.size());
		RSort_Parallel(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size(), s.data());

		mortonTree[0].range = { 0, (int)pointsCopy.size() - 1 };
		mortonTree[0].particle = false;
	}//MakeRootMortonTree()


	//Функция вычисления длины общей части префиксов двух(а значит - и диапазона) частиц
	int MortonTree::Delta(int i, int j) const
	{
		if ((j < 0) || (j > (int)pointsCopy.size() - 1))
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
	}//Delta(...)


	//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
	int MortonTree::PrefixLength(int cell) const
	{
		const auto& range = mortonTree[cell].range;
		return std::min(Delta(range[0], range[1]), 2 * codeLength);
	}//PrefixLength(...)

	
	//Функция вычисления общего префикса двух частиц
	std::pair<unsigned int, int> MortonTree::Prefix(int cell) const
	{
		int length = PrefixLength(cell);
		unsigned int el = mortonCodes[mortonTree[cell].range[0]].key;
		return { el >> (2 * codeLength - length), length };
	}//Prefix(...)


	//Функция вычисления геометрических параметров внутренней ячейки дерева
	void MortonTree::SetCellGeometry(int cell)
	{
		int prLength;
		unsigned int pr;
		std::tie<unsigned int, int>(pr, prLength) = Prefix(cell);
		prLength -= PrefixLength(0);

		Point2D sz = { 1.0 / (double)(1 << ceilhalf(prLength)), 1.0 / (1 << (prLength / 2)) };

		Point2D pos = 0.5 * sz;

		for (int i = 0; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[0] += 1.0 / (1 << (1 + i / 2));

		for (int i = 1; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[1] += 1.0 / (1 << (1 + i / 2));

		mortonTree[cell].center = pos * (1.0 / iQuadSideVar) + lowLeft;
		mortonTree[cell].size = sz * (1.0 / iQuadSideVar);
		mortonTree[cell].level = prLength;

		mortonTree[cell].leaf = (prLength >= maxTreeLevel);
	}//SetCellGeometry(...)


	//Функция построения i-й внутренней ячейки дерева
	void MortonTree::BuildInternalTreeCell(int i)
	{
		const int n = (int)pointsCopy.size();
		int d = sign(Delta(i, i + 1) - Delta(i, i - 1));
		int delta_min = Delta(i, i - d);
		
		int Lmax = 2;
		while (Delta(i, i + Lmax * d) > delta_min)
			Lmax *= 2;

		int L = 0;
		for (int t = (Lmax >> 1); t >= 1; t >>= 1)		
			if (Delta(i, i + (L + t) * d) > delta_min)
				L += t;
		
		int j = i + L * d;
		int delta_node = Delta(i, j);
		int s = 0;
		for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
		{
			int dl = Delta(i, i + (s + t) * d);
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
	}//BuildInternalTreeCell(...)


	//Построение внутренних ячеек дерева
	void MortonTree::BuildMortonInternalTree()
	{
#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < (int)pointsCopy.size() - 1; ++q)
			BuildInternalTreeCell(q);
		
		mortonTree[pointsCopy.size() - 1].leaf = false;
		mortonTree[pointsCopy.size() - 1].particle = false;

#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < (int)pointsCopy.size() - 1; ++q)
			SetCellGeometry(q);
	}//BuildMortonInternalTree()


	//Построение верхушек дерева --- отдельных частиц
	void MortonTree::BuildMortonParticlesTree()
	{
#pragma omp parallel for		
		for (int q = 0; q < (int)pointsCopy.size(); ++q)
		{
			auto& pt = mortonCodes[q];

			auto& cell = mortonTree[pointsCopy.size() + q];
			cell.center = pt.r;

			cell.leaf = true;
			cell.size = { 0.0, 0.0 };
		}//for q
	}//BuildMortonParticlesTree(...)


	//Заполнение списка нижних вершин: 
	//рекурсивный алгоритм, быстро работает в последовательном варианте		
	void MortonTree::FillMortonLowCells(int cell)
	{
		if (mortonTree[cell].leaf)
			mortonLowCells.push_back(cell);
		else
		{
			FillMortonLowCells(mortonTree[cell].child[0]);
			FillMortonLowCells(mortonTree[cell].child[1]);
		}//else
	}//FillMortonLowCells(...)


	//Заполнение списка нижних вершин: 
	//нерекурсивный алгоритм, хорошо распараллеливается, но неэффективен в последовательном варианте
	void MortonTree::FillMortonLowCellsA()
	{
#pragma omp parallel
		{
			std::vector<int> good_matches_private;
#pragma omp for nowait
			for (int cell = 0; cell < (int)mortonTree.size(); ++cell)			
				if (!mortonTree[mortonTree[cell].parent].leaf && mortonTree[cell].leaf)				
					good_matches_private.push_back(cell);				
#pragma omp critical
			mortonLowCells.insert(mortonLowCells.end(), good_matches_private.begin(), good_matches_private.end());
		}//parallel
	}//FillMortonLowCellsA()


	/// Вычисление параметров дерева (циркуляций и высших моментов)
	void MortonTree::CalculateMortonTreeParams(int cell, int omplevel)
	{
		auto& cl = mortonTree[cell];
		if (!cl.particle)
		{
#pragma omp parallel /*default(none)*/ shared(cl, omplevel) num_threads(2) if (omplevel < prm.maxLevelOmp + 1)
			{				
				std::vector<Point2D> h(prm.order);
#pragma omp for 
				for (int s = 0; s < 2; ++s)
					CalculateMortonTreeParams(cl.child[s], omplevel + 1);

#pragma omp single
				{					
					for (auto& m : cl.mom)
						m.toZero();
										
					for (int s = 0; s < 2; ++s)//цикл по двум потомкам
					{
						auto& chld = mortonTree[cl.child[s]];
						const Point2D& g = chld.mom[0];
						cl.mom[0] += g;

						
						h[0] = (chld.center) - (cl.center);
						for (int i = 1; i < prm.order; ++i)
							h[i] = multz(h[i - 1], h[0]);
						
						for (int p = 1; p <= prm.order; ++p) 
						{
							Point2D shiftMom = chld.mom[p];
							for (int k = 1; k <= p; ++k)
							{
								shiftMom += getBinom(p, k) *  multz(chld.mom[p - k], h[k - 1]);
								ADDOP(2);
							}//for k
							cl.mom[p] += shiftMom;
						}//for m
						
					}//for s
				}//single				
			}//parallel
		}//if(!particle)
		else
		{
			//расчет моментов ячейки нижнего уровня для вихрей или для панелей
			if (pans) //если дерево из панелей
			{
				const auto & panel = pointsCopy[mortonCodes[cl.range[0]].originNumber];
				for (int p = 0; p <= prm.order; p += 2)
					cl.mom[p] = panel.Mpan[p] * panel.g();
				ADDOP(prm.order/2);
#ifdef linScheme
				for (int p = 1; p <= prm.order; p += 2)
					cl.mom[p] = panel.Mpan[p] * panel.gamLin;
#ifdef asympScheme
				if (panel.MpanAs.size() > 0)
				{					
					for (int p = 1; p <= prm.order; ++p)
						cl.mom[p] += panel.MpanAs[p] * panel.gamLin;  
					ADDOP(prm.order/2);
				}
#endif //asympScheme
#else
				for (int p = 1; p <= prm.order; p += 2)
					cl.mom[p].toZero();//Для константной схемы
#endif // linScheme
				
			}
			else //если дерево из частиц
			{
				//C учетом того, что дерево строится до частиц --- только монопольный момент отличен от нуля
				cl.mom[0][0] = mortonCodes[cl.range[0]].g;
				cl.mom[0][1] = 0.0;

				for (int p = 1; p <= prm.order; ++p)
					cl.mom[p].toZero();
			}//else
		}//else
	}//CalculateMortonTreeParams(...)		

	//Расчет коэффициентов разложения в ряд Тейлора внутри ячейки нижнего уровня
	void MortonTree::CalcLocalCoeffToLowLevel(int lowCell, std::unique_ptr<MortonTree>& treeInf, int fromWho, bool calcCloseTrees)
	{		
		auto& lt = mortonTree[lowCell];
		//if (calcCloseTrees)
		//	lt.closeCells.resize(0);

		if (treeInf.get() != this || lowCell != fromWho)
		{
			auto& wh = treeInf->mortonTree[fromWho];

			double h, h0;
			h = wh.size[0] + wh.size[1];
			h0 = lt.size[0] + lt.size[1];

			Point2D& r0 = lt.center;
			Point2D& r1 = wh.center;

			Point2D rr = r0 - r1;

			double crit2 = rr.length2();
			ADDOP(3);
			ADDOP(3);
			// если выполнен критерий дальности => считаем коэффициенты
			if ((crit2 >= sqr((h0 + h + 2.0 * prm.eps) / prm.theta)))
			{
				Point2D  rr = r0 - r1;
				rr *= 1.0 / rr.length2();
				ADDOP(3);
				
				Point2D varTheta = rr;
				for (int diag = 0; diag < prm.order; ++diag)
				{
					for (int q = 0; q <= diag; ++q)
					{
						lt.E[q] += ((q % 2 ? -1.0 : 1.0) * iFact[diag - q]) * multzA(varTheta, wh.mom[diag - q]);
						ADDOP(3);
					}
					varTheta = (diag + 1) * multz(varTheta, rr);
					ADDOP(2);
				}
				lt.E[0] += iFact[prm.order] * multzA(varTheta, wh.mom[prm.order]);
				ADDOP(2);

			}//if crit2
			else // если не выполнен критерий, то рекурсия 
			{
				if (!wh.leaf)
				{
					CalcLocalCoeffToLowLevel(lowCell, treeInf, wh.child[0], calcCloseTrees);
					CalcLocalCoeffToLowLevel(lowCell, treeInf, wh.child[1], calcCloseTrees);
				}
				else if (calcCloseTrees)
				{					
					//std::cout << fromWho << std::endl;
					if (treeInf->index < 0)
						lt.closeCells.push_back(fromWho);
					else
						lt.closeCellsPfl[treeInf->index].push_back(fromWho);
				}
			}//else if
		}//if (lowTree != this)
		else if (calcCloseTrees)
		{
			if (treeInf->index < 0)
				lt.closeCells.push_back(lowCell); //себя тоже добавляем в ближнюю зону 
			else
				lt.closeCellsPfl[treeInf->index].push_back(lowCell);
		}
	}//CalcLocalCoeffToLowLevel(...)



	//Расчет влияния от ближних ячеек по формуле Био - Савара
	void MortonTree::CalcVeloBiotSavart(int lowCell, std::unique_ptr<MortonTree>& treeInf)
	{
		auto& lc = mortonTree[lowCell];
		for (size_t k = 0; k < lc.closeCells.size(); ++k)
		{
			//Локальные переменные для цикла
			Point2D velI, velLinI;
			double dst2eps;

			for (int i = lc.range[0]; i <= lc.range[1]; ++i)
			{
				velI.toZero();
				velLinI.toZero();

				const Point2D& posI = mortonCodes[i].r; 

				auto& rg = treeInf->mortonTree[lc.closeCells[k]].range;
				for (int j = rg[0]; j <= rg[1]; ++j)
				{
					auto& pt = treeInf->mortonCodes[j];

					if (!pans) //В контрольном дереве - точки
					{
						if (!treeInf->pans)// если точки влияют на точки
						{
							const Point2D& posJ = pt.r;
							const double& gamJ = pt.g;
							dst2eps = std::max((posI - posJ).length2(), prm.eps2);
							velI += (gamJ / dst2eps) * (posI - posJ).kcross();
							ADDOP(5);
						}
						else // панели влияют на точки
						{
							const auto& pnl = treeInf->pointsCopy[pt.originNumber];
							Point2D u0 = pnl.tau;
							Point2D pp = posI - pnl.panEnd;
							Point2D s = posI - pnl.panBegin;
							double alpha = atan2(pp ^ s, pp & s);
							double lambda = 0.5 * log((s & s) / (pp & pp));
							ADDOP(11);

							Point2D va = pp + s;
							Point2D vb = pnl.tau;
							Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
							Point2D u1 = (0.5 / pnl.len) * omega;

							velI += (pnl.g() / pnl.len) * (-alpha * (u0.kcross()) + lambda * u0).kcross();
							ADDOP(18);
#ifdef linScheme
							velI += (pnl.gamLin / pnl.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnl.tau).kcross();
							ADDOP(7);
#ifdef asympScheme						
							for (int t = 0; t < prm.nAngPoints; ++t)
							{
								if (pt.originNumber == prm.KK[t])
								{
									const auto& pnl = treeInf->pointsCopy[prm.KK[t]];
									Point2D u0 = pnl.tau;
									Point2D pp = posI - pnl.panEnd;
									Point2D s = posI - pnl.panBegin;
									double alpha = atan2(pp ^ s, pp & s);
									double lambda = 0.5 * log((s & s) / (pp & pp));
									ADDOP(12);

									Point2D va = pp + s;
									Point2D vb = pnl.tau;
									Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
									Point2D u1 = (0.5 / pnl.len) * omega;
									ADDOP(10);

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
									ADDOP(35);
								}
								if (pt.originNumber == prm.KKm[t])
								{
									//auto& pnl = pointsCopyPan[KKm[t]];
									//const Point2D& posI = { 1, 1 };
									const auto& pnl = treeInf->pointsCopy[prm.KKm[t]];
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
									ADDOP(35);
								}
							}
#endif//asympScheme
#endif//linScheme	
						}
					}//if (!pans)
					else//В контрольном дереве - панели
					{
						const Point2D& posJ = pt.r;
						const auto& pnlI = pointsCopy[mortonCodes[i].originNumber];
					
						
						Point2D u0 = pnlI.tau;
						Point2D pp = posJ - pnlI.panEnd;
						Point2D s = posJ - pnlI.panBegin;
						double alpha = atan2(pp ^ s, pp & s);
						double lambda = 0.5 * log((s & s) / (pp & pp));

						Point2D va = pp + s;
						Point2D vb = pnlI.tau;
						Point2D omega = (va & vb) * vb + (va ^ vb) * (vb.kcross());
						Point2D u1 = (0.5 / pnlI.len) * omega;

						velI -= (pt.g / pnlI.len) * (-alpha * (u0.kcross()) + lambda * u0).kcross();
						ADDOP(30);
#ifdef linScheme
						velLinI -= (pt.g / pnlI.len) * (-alpha * (u1.kcross()) + lambda * u1 - pnlI.tau).kcross();
						ADDOP(7);
#endif
					}//else if (!pans)

				}//for j

				pointsCopy[mortonCodes[i].originNumber].veloCopy += velI;
#ifdef linScheme
				pointsCopy[mortonCodes[i].originNumber].veloCopyLin += velLinI;
#endif
			}//for i 
		}//for k
	}//CalcVeloBiotSavart(...)

	//Расчет интегрального влияния от распределения завихренности панелей ближних ячеек (для решения СЛАУ)
	void MortonTree::CalcInfluenceFromPanels(int lowCell, std::unique_ptr<MortonTree>& treeInf)
	{
		if ((treeInf->index == 0))
		{
			auto& lc = mortonTree[lowCell];
			for (int i = lc.range[0]; i <= lc.range[1]; ++i)
			{
				PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];
				itDop.i00.resize(prm.airfoilFile.size());
#ifdef linScheme
				itDop.i01.resize(prm.airfoilFile.size());
				itDop.i10.resize(prm.airfoilFile.size());
				itDop.i11.resize(prm.airfoilFile.size());
#endif
			}
		}

		numvector<double, 3> alpha, lambda;

		//auxillary vectors
		Point2D p1, s1, p2, s2, di, dj, i00, i01, i10, i11;
		numvector<Point2D, 3> v00, v11;
		numvector<Point2D, 2> v01, v10;
		
		numvector<Point2D, 2> H01, H11;

		int n = (int)pointsCopy.size();
		auto& lc = mortonTree[lowCell];

		//std::cout << "lc.closeCells.size() = " << lc.closeCells.size() << std::endl;

		for (size_t k = 0; k < lc.closeCellsPfl[treeInf->index].size(); ++k)
		{
			//Локальные переменные для цикла
			Point2D velI, velIlin;

			const auto& rg = treeInf->mortonTree[lc.closeCellsPfl[treeInf->index][k]].range;

			//#pragma omp parallel for schedule(dynamic,10) ordered
			//Обязательно должно быть закрыто распараллеливание!
			for (int i = lc.range[0]; i <= lc.range[1]; ++i)
			{
				PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];

				velI.toZero();
				velIlin.toZero();

				const Point2D& posI = itDop.r();
				di = itDop.panEnd - itDop.panBegin;
				double dilen = di.length();
				double ileni = 1.0 / dilen;

				const Point2D& taui = ileni * di;
				ADDOP(6);

				const double cftReserve = 0.01;
				if (k == 0)
				{					
					itDop.i00[treeInf->index].resize(0);
					itDop.i00[treeInf->index].reserve(int(cftReserve * n));

#ifdef linScheme
					itDop.i01[treeInf->index].resize(0);
					itDop.i10[treeInf->index].resize(0);
					itDop.i11[treeInf->index].resize(0);
					itDop.i01[treeInf->index].reserve(int(cftReserve * n));
					itDop.i10[treeInf->index].reserve(int(cftReserve * n));
					itDop.i11[treeInf->index].reserve(int(cftReserve * n));
#endif
				}

				for (int j = rg[0]; j <= rg[1]; ++j)
				{
					auto& pt = treeInf->pointsCopy[treeInf->mortonCodes[j].originNumber]; //(mortonCodes[j].originNumber)-я панель
					const Point2D& posJ = pt.r();
					const double& gamJ = pt.g();
#ifdef linScheme
					const double& gamJlin = pt.gamLin;
#endif
					if (posI != posJ)
					{
						dj = pt.panEnd - pt.panBegin;
						double djlen = dj.length();
						double ilenj = 1.0 / djlen;
						const Point2D& tauj = ilenj * dj;

						p1 = itDop.panEnd - pt.panEnd;
						s1 = itDop.panEnd - pt.panBegin;
						p2 = itDop.panBegin - pt.panEnd;
						s2 = itDop.panBegin - pt.panBegin;

						alpha = { \
									isAfter(pt, itDop) ? 0.0 : Alpha(s2, s1), \
									Alpha(s2, p1), \
									isAfter(itDop, pt) ? 0.0 : Alpha(p1, p2) \
						};

						lambda = { \
									isAfter(pt, itDop) ? 0.0 : Lambda(s2, s1), \
									Lambda(s2, p1), \
									isAfter(itDop, pt) ? 0.0 : Lambda(p1, p2) \
						};

						v00 = { \
								Omega(s1, taui, tauj), \
								- Omega(di, taui, tauj), \
								Omega(p2, taui, tauj) \
						};

						i00 = ileni * ilenj * (alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2] \
							+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]).kcross());
						ADDOP(21);
#ifdef linScheme
						double s1len2 = s1.length2();
						ADDOP(2);

						v01 = {
								0.5 / (djlen) * (((p1 + s1) & tauj) * Omega(s1, taui, tauj) - s1len2 * taui),
								-0.5 * dilen / djlen * Omega(s1 + p2, tauj, tauj)
						};

						i01 = ileni * ilenj * ((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]\
							+ (((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * dilen * tauj).kcross());
						ADDOP(27);

						v10 = {
								-0.5 / dilen * (((s1 + s2) & taui) * Omega(s1, taui, tauj) - s1len2 * tauj),
								0.5 * djlen / dilen * Omega(s1 + p2, taui, taui)
						};

						i10 = ileni * ilenj * ((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1] \
							+ (((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * djlen * taui).kcross());
						ADDOP(27);


						v11 = {
							1.0 / (12.0 * dilen * djlen) * \
							(2.0 * (s1 & Omega(s1 - 3.0 * p2, taui, tauj)) * Omega(s1, taui, tauj) - \
									s1len2 * (s1 - 3.0 * p2)) - 0.25 * Omega(s1, taui, tauj),
							-dilen / (12.0 * djlen) * Omega(di, tauj, tauj),
							-djlen / (12.0 * dilen) * Omega(dj, taui, taui)
						};


						i11 = ileni * ilenj * ((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]\
							+ ((lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
								+ 1.0 / 12.0 * (djlen * taui + dilen * tauj - 2.0 * Omega(s1, taui, tauj))).kcross());
						ADDOP(50);

#ifdef asympScheme
						for (int t = 0; t < prm.nAngPoints; t++) {
							if (mortonCodes[j].originNumber == prm.KK[t])
							{

								H01[0].toZero();
								H01[1].toZero();

								//auto H01old = StoConstFrom1({ pt.panBegin, pt.panEnd }, { itDop.panBegin, itDop.panEnd }, t, prm);
								if (mortonCodes[i].originNumber == prm.KK[t] + 1)
									H01 = StoConstFrom1_new1({ pt.panBegin, pt.panEnd }, { itDop.panBegin, itDop.panEnd }, t, prm);
								else
								if (mortonCodes[i].originNumber == prm.KKm[t])
									H01 = StoConstFrom1_new2({ pt.panBegin, pt.panEnd }, { itDop.panBegin, itDop.panEnd }, t, prm);
								else
									H01 = StoConstFrom1_new({ pt.panBegin, pt.panEnd }, { itDop.panBegin, itDop.panEnd }, t, prm);

								i01 = /*taui &*/ (ilenj * ileni * DPI * (H01[0].kcross())) - (1.0 / (1.0 - prm.mu[t]) * i00);
								i11 = /*taui &*/ (ilenj * ileni * DPI * ((H01[1]-0.5 * H01[0]).kcross())) - (1.0 / (1.0 - prm.mu[t]) * i10);
								ADDOP(16);
							}

							if (mortonCodes[j].originNumber == prm.KKm[t])
							{
								H11[0].toZero();
								H11[1].toZero();


								if (mortonCodes[i].originNumber == prm.KK[t])
										H11 = StoConstFrom1_new2({ pt.panEnd, pt.panBegin }, { itDop.panEnd, itDop.panBegin }, t, prm);
								else
									if (mortonCodes[i].originNumber == prm.KKm[t] - 1) 
										H11 = StoConstFrom1_new1({ pt.panEnd, pt.panBegin }, { itDop.panEnd, itDop.panBegin }, t, prm);
									else
										H11 = StoConstFrom1_new({ pt.panEnd, pt.panBegin }, { itDop.panEnd, itDop.panBegin }, t, prm);

								i01 = ileni * ilenj * DPI * H11[0].kcross() - (1.0 / (1.0 - prm.mu[t]) * i00);
								i11 = -(ileni * ilenj * DPI * (H11[1] - 0.5 * H11[0]).kcross() + (1.0 / (1.0 - prm.mu[t]) * i10));
								ADDOP(16);
							}
						}
#endif

#endif
						if (treeInf.get() == this)
						{
							if (isAfter(itDop, pt))
							{
								itDop.a = i00 * IDPI;
								ADDOP(2);

#ifdef linScheme
								itDop.a1 = i11 * IDPI;
								ADDOP(2);
#endif
							}
							if (isAfter(pt, itDop))
							{
								itDop.c = i00 * IDPI;
								ADDOP(2);

#ifdef linScheme
								itDop.c1 = i11 * IDPI;
								ADDOP(2);
#endif
							}
						}
//#pragma omp ordered
{
						SizeCheck(itDop.i00[treeInf->index]);
						itDop.i00[treeInf->index].push_back(i00);

						//if (itDop.i00.size() > 5000)
						//	std::cout << "itDop.i00.size = " << itDop.i00.size() << std::endl;

						velI += gamJ * i00;
						ADDOP(2);

#ifdef linScheme					

						SizeCheck(itDop.i01[treeInf->index]);
						itDop.i01[treeInf->index].push_back(i01);

						SizeCheck(itDop.i10[treeInf->index]);
						itDop.i10[treeInf->index].push_back(i10);

						SizeCheck(itDop.i11[treeInf->index]);
						itDop.i11[treeInf->index].push_back(i11);

						velI += gamJlin * i01;
						velIlin += gamJ * i10 + gamJlin * i11;
						ADDOP(6);

#endif
						}
					}//if (posI != posJ)
				}//for j

				itDop.veloCopy += velI;
#ifdef linScheme
				itDop.veloCopyLin += velIlin;
#endif
			}//for i 
		}//for k
	}//CalcInfluenceFromPanels(...)


	// Обновление влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
	void MortonTree::UpdateInfluence(int lowCell, std::unique_ptr<MortonTree>& treeInf)
	{
		auto& lc = mortonTree[lowCell];

		std::vector<int> s(lc.range[1] - lc.range[0] + 1, 0);

		for (size_t k = 0; k < lc.closeCellsPfl[treeInf->index].size(); ++k)
		{
			//Локальные переменные для цикла
			Point2D velI, velIlin;			
			auto& rg = treeInf->mortonTree[lc.closeCellsPfl[treeInf->index][k]].range;
			int inum = 0;

			//if (lc.range[1] - lc.range[0] > 0)
			//	std::cout << "!!!" << std::endl;

			for (int i = lc.range[0]; i <= lc.range[1]; ++i, ++inum)
			{
				PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];

				velI.toZero();
				velIlin.toZero();

				const Point2D& posI = itDop.r();
				
				//if (rg[1]-rg[0] > 0)
				//	std::cout << "???" << std::endl;

				for (int j = rg[0]; j <= rg[1]; ++j)
				{
					const PointsCopy& itDop2 = treeInf->pointsCopy[mortonCodes[j].originNumber];
					const Point2D& posJ = itDop2.r();
					const double& gamJ = itDop2.g();
#ifdef linScheme	
					const double& gamJlin = itDop2.gamLin;
#endif

					if (posI != posJ)
					{
						velI += gamJ * (itDop.i00[treeInf->index][s[inum]]);
						ADDOP(2);

#ifdef linScheme					
						velI += gamJlin * (itDop.i01[treeInf->index][s[inum]]);
						velIlin += gamJ * (itDop.i10[treeInf->index][s[inum]]) + gamJlin * (itDop.i11[treeInf->index][s[inum]]);
						ADDOP(6);

#endif
						++s[inum];
					}
				}//for j
			itDop.veloCopy += velI;
#ifdef linScheme	
			itDop.veloCopyLin += velIlin;
#endif
			}//for i
		}//for k
	}//UpdateInfluence(...)


	//Расчет точечного влияния от распределения завихренности панелей ближних ячеек (для решения СЛАУ)
	void MortonTree::CalcInfluenceFromPanelsToPoints(int lowCell)
	{

	}//CalcInfluenceFromPanelsToPoints(...)

	//Вычисляет влияния части подряд идущих вихрей из вихревого следа на прямолинейную панель для правой части
	void GetInfluenceFromVorticesToPanel(PointsCopy panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs)
	{
		double& velI = wakeRhs[0];

		const Point2D& posI0 = panel.panBegin;
		const Point2D& posI1 = panel.panEnd;

#ifdef linScheme
		double& velILin = wakeRhs[1];
		const Point2D& taui = panel.tau;
		Point2D di = posI1 - posI0;
#endif // linScheme

		for (size_t it = 0; it != count; ++it)
		{
			const Vortex2D& vt = ptr[it];
			const Point2D& posJ = vt.r();
			const double& gamJ = vt.g();

			Point2D s = posJ - posI0;
			Point2D p = posJ - posI1;

			double alpha = Alpha(p, s);
			velI -= gamJ * alpha;
			ADDOP(2);

#ifdef linScheme
			Point2D u1 = 0.5 / di.length() * Omega(p + s, taui, taui);
			double lambda = Lambda(p, s);
			velILin -= gamJ * (alpha * (u1 & taui) + lambda * (u1 & (-taui.kcross())));
			ADDOP(16);

#endif // linScheme
		}
	}//GetInfluenceFromVorticesToPanel(...)





	//Расчет влияния от следа на панели профиля
	void MortonTree::CalcInfluenceFromPointsToPanel(int lowCell, std::unique_ptr<MortonTree>& treeInf)
	{
		auto& lc = mortonTree[lowCell];
		std::vector<Vortex2D> vtxClose;

		//собираем массив ближних вихрей для ячейки lc
		for (size_t k = 0; k < lc.closeCells.size(); ++k)
		{
			auto& rg = treeInf->mortonTree[lc.closeCells[k]].range;
			for (int j = rg[0]; j <= rg[1]; ++j)
				vtxClose.push_back(Vortex2D(treeInf->mortonCodes[j].r, treeInf->mortonCodes[j].g));
		}

		size_t shDim = 1;
#ifdef linScheme
		shDim = 2;
#endif // linScheme

		std::vector<double> velI(shDim, 0.0);

		for (int i = lc.range[0]; i <= lc.range[1]; ++i)
		{
			auto& panel = pointsCopy[mortonCodes[i].originNumber];
			GetInfluenceFromVorticesToPanel(panel, vtxClose.data(), vtxClose.size(), velI);

			for (size_t j = 0; j < shDim; ++j)
			{
				velI[j] *= IDPI / panel.len;
				ADDOP(3);
			}
			panel.velTau = velI[0];
#ifdef linScheme
			panel.velTauLin = velI[1];
#endif // linScheme
		}//for i
	}//CalcInfluenceFromPointsToPanel(...)



	// Расчет влияния от дальней зоны при помощи Тейлоровского разложения
	void MortonTree::CalcVeloTaylorExpansion(int lowCell)
	{
		auto& lc = mortonTree[lowCell];
		for (int i = lc.range[0]; i <= lc.range[1]; ++i)
		{
			Point2D deltaPos = mortonCodes[i].r - lc.center;
			Point2D v = lc.E[0];

#ifndef infToPanels
			Point2D hh = deltaPos;
			for (int k = 1; k < prm.order - 1; ++k) 
			{
				v += iFact[k] * multzA(lc.E[k], hh);
				hh = multz(hh, deltaPos);
				ADDOP(2);
			}//for k
			v += iFact[prm.order-1] * multzA(lc.E[prm.order-1], hh);
			ADDOP(2);
#else
			PointsCopy& itDop = pointsCopy[mortonCodes[i].originNumber];

			Point2D rPan = itDop.panEnd - itDop.panBegin;
			//std::cout << "rpan = " << rPan << std::endl;
			Point2D dPos2 = 2.0 * deltaPos;

			Point2D kp, km;

			Point2D mulP = kp = dPos2 + rPan;
			Point2D mulM = km = dPos2 - rPan;

			Point2D taudL = (0.5 / itDop.len) * itDop.tau;
			ADDOP(5);


#ifdef linScheme
			Point2D vL;
			vL.toZero();
			Point2D taudLc = taudL;
#endif
			for (int k = 1; k < prm.order; ++k)
			{
				mulP = multz(mulP, kp);
				mulM = multz(mulM, km);
				taudL /= 2.0;
				v += iFact[k + 1] * multz(lc.E[k], multzA(taudL, mulP - mulM));	
				ADDOP(2);

#ifdef linScheme
				vL += (iFact[k + 1] / (k + 2)) * multz(lc.E[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan - dPos2) + multz(mulM, (k + 1) * rPan + dPos2)));
				ADDOP(3);

#endif
			}
#endif // infToPanels
			pointsCopy[mortonCodes[i].originNumber].veloCopy += Point2D({ -v[1], v[0] });
	
#ifdef infToPanels
#ifdef linScheme
			pointsCopy[mortonCodes[i].originNumber].veloCopyLin += Point2D({ -vL[1], vL[0] });
#endif
#endif
		}//for i
	}//CalcVeloTaylorExpansion(...)
	


	void MortonTree::AddVelo(int lowCell)
	{
		auto& lc = mortonTree[lowCell];
		for (int j = lc.range[0]; j <= lc.range[1]; ++j)
		{
			auto& pnt = pointsCopy[mortonCodes[j].originNumber];
			pnt.velTau += pnt.veloCopy & pnt.tau;
			pnt.velTauLin += pnt.veloCopyLin & pnt.tau;
			ADDOP(4);
		}
	}//AddVelo(...)



	//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
	unsigned int MortonTree::ExpandBits(unsigned int v)
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
	unsigned int MortonTree::Morton2D(const Point2D& r)
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
		const unsigned int& xx = ExpandBits((unsigned int)(rscale[0]));
		const unsigned int& yy = ExpandBits((unsigned int)(rscale[1]));
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
		if (n == 0)
			return;
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

	void MortonTree::getStatistics(int& numParticles, int& treshold, std::pair<int, int>& partLevel, std::pair<int, int>& leafsLevel, int& lowLevelCells) const
	{
		numParticles = (int)pointsCopy.size();
		treshold = maxTreeLevel;

		int leafsLevelMin = (int)1e+9;
		int leafsLevelMax = 0;
		for (int q = 0; q < numParticles - 1; ++q)
		{
			if (!mortonTree[mortonTree[q].parent].leaf)
			{
				if (mortonTree[q].leaf)
				{
					int level = mortonTree[q].level;

					if (level > leafsLevelMax)
						leafsLevelMax = level;

					if (level < leafsLevelMin)
						leafsLevelMin = level;

					continue;
				}

				if ((mortonTree[mortonTree[q].child[0]].particle) || (mortonTree[mortonTree[q].child[1]].particle))
				{
					int level = mortonTree[q].level + 1;
					if (level > leafsLevelMax)
						leafsLevelMax = level;

					if (level < leafsLevelMin)
						leafsLevelMin = level;

					continue;
				}
			}
		}
		leafsLevel = { leafsLevelMin, leafsLevelMax };


		int partLevelMin = (int)1e+9;
		int partLevelMax = 0;
		for (int q = 0; q < numParticles - 1; ++q)
		{
			if ((mortonTree[mortonTree[q].child[0]].particle) || (mortonTree[mortonTree[q].child[1]].particle))
			{
				int level = mortonTree[q].level;
				if (level > partLevelMax)
					partLevelMax = level;

				if (level < partLevelMin)
					partLevelMin = level;
			}
		}
		partLevel = { partLevelMin + 1, partLevelMax + 1 };
		lowLevelCells = (int)mortonLowCells.size();
	}


}//namespace BH

