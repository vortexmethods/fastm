/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.0    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2021/08/05     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
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
\version 1.0
\date 05 августа 2021 г.
*/

#include <algorithm>

#include <omp.h>

#include "Tree.h"

namespace BH
{

	const double PI = 3.1415926535897932384626;
	const double IPI = 1.0 / PI;
	const double IDPI = 0.5 / PI;

	extern long long op;


	///Построение корня дерева на основе заданных вихрей 
	void Tree::makeTree(std::vector<PointsCopy>& points)
	{
		auto minmaxY = std::minmax_element(points.begin(), points.end(), [](const PointsCopy& a, const PointsCopy& b)
			{
				return a.r()[1] < b.r()[1];
			});

		iBottom = (minmaxY.first)->r()[1];
		iTop = (minmaxY.second)->r()[1];

		auto minmaxX = std::minmax_element(points.begin(), points.end(), [](const PointsCopy& a, const PointsCopy& b)
			{
				return a.r()[0] < b.r()[0];
			});

		iLeft = (minmaxX.first)->r()[0];
		iRight = (minmaxX.second)->r()[0];

		index.resize(points.size());
		for (int i = 0; i < points.size(); ++i)
			index[i] = i;

		posCentre = { 0.5 * (iLeft + iRight), 0.5 * (iBottom + iTop) };

		level = 0;
		p = 0;

		treePtr.resize((1 << (NumOfLevels + 1)) - 1);

		for (auto& p : treePtr)
			p.reset(new Tree(points));
	}//makeTree(...)


	/// Построение структуры дерева 
	void Tree::CreateTree()
	{
		auto lamx = [this](int i, int j)
		{
			return points[i].r()[0] < points[j].r()[0];
		};

		auto lamy = [this](int i, int j)
		{
			return points[i].r()[1] < points[j].r()[1];
		};


#pragma omp parallel default(none) shared(std::cout, lamx, lamy) num_threads(1) if(level <= maxLevelOmp)
		{
#pragma omp single
			{
				double w = iRight - iLeft;
				double h = iTop - iBottom;

				int k2 = 1 << (level + 1);

				mChildren[0] = std::move(BigTree->treePtr[k2 + 2 * p - 1]);
				mChildren[1] = std::move(BigTree->treePtr[k2 + 2 * p]);

				if (w > h)
				{
					double midX = 0.5 * (iLeft + iRight);

					mChildren[0]->iLeft = iLeft;
					mChildren[0]->iRight = midX;
					mChildren[0]->iBottom = iBottom;
					mChildren[0]->iTop = iTop;

					mChildren[1]->iLeft = midX;
					mChildren[1]->iRight = iRight;
					mChildren[1]->iBottom = iBottom;
					mChildren[1]->iTop = iTop;

					auto itp = std::partition(index.begin(), index.end(), [this, midX](int a) {return points[a].r()[0] < midX; });

					mChildren[0]->index.insert(mChildren[0]->index.end(), index.begin(), itp);
					mChildren[1]->index.insert(mChildren[1]->index.end(), itp, index.end());
				} //if (w > h)
				else
				{
					double midY = 0.5 * (iBottom + iTop);

					mChildren[0]->iLeft = iLeft;
					mChildren[0]->iRight = iRight;
					mChildren[0]->iBottom = iBottom;
					mChildren[0]->iTop = midY;

					mChildren[1]->iLeft = iLeft;
					mChildren[1]->iRight = iRight;
					mChildren[1]->iBottom = midY;
					mChildren[1]->iTop = iTop;

					auto itp = std::partition(index.begin(), index.end(), [this, midY](int a) {return points[a].r()[1] < midY; });

					mChildren[0]->index.insert(mChildren[0]->index.end(), index.begin(), itp);
					mChildren[1]->index.insert(mChildren[1]->index.end(), itp, index.end());
				}//else
			}

#pragma omp for 
			for (int s = 0; s < 2; ++s)
			{
				auto& mC = *mChildren[s];

				//обрезаем по вихрям
				auto minmaxX = std::minmax_element(mC.index.begin(), mC.index.end(), lamx);

				mC.iLeft = points[*(minmaxX.first)].r()[0];
				mC.iRight = points[*(minmaxX.second)].r()[0];

				auto minmaxY = std::minmax_element(mC.index.begin(), mC.index.end(), lamy);

				mC.iBottom = points[*(minmaxY.first)].r()[1];
				mC.iTop = points[*(minmaxY.second)].r()[1];

				mC.level = level + 1;
				mC.p = 2 * p + s;

				mC.BigTree = BigTree;

				mC.posCentre = { 0.5 * (mC.iLeft + mC.iRight), 0.5 * (mC.iBottom + mC.iTop) };


				if ((mC.level < NumOfLevels) && (mC.index.size() > minNumOfVort) && \
					(std::max(fabs(mC.iRight - mC.iLeft), fabs(mC.iTop - mC.iBottom)) > minDist))
				{
					mC.CreateTree();
				}
				else
				{
					mC.lowLevel = true;
					mC.closeTrees.reserve(100);
				}
			}
		}
	}//CreateTree(...)


	void Tree::PushbackLowTrees()
	{
		if (!lowLevel)
			for (int s = 0; s < 2; ++s)
				mChildren[s]->PushbackLowTrees();
		else
			BigTree->lowTrees.push_back(this);

	}//PushbackLowTrees()


	/// Вычисление параметров дерева (циркуляций и центров завихренности)
	void Tree::CalculateTreeParams()
	{
		if (!lowLevel)
		{
#pragma omp parallel default(none) num_threads(2) if(level <= maxLevelOmp)
			{
#pragma omp for 
				for (int s = 0; s < 2; ++s)
					mChildren[s]->CalculateTreeParams();

#pragma omp single
				{
					mon = 0.0;
					dip.toZero();
					qua.toZero();
					oct.toZero();
					hex.toZero();

					Point2D h;
					for (int s = 0; s < 2; ++s)
					{
						const double& g = mChildren[s]->mon;
						mon += g;

						h = (mChildren[s]->posCentre) - (this->posCentre);

						const Point2D& d = mChildren[s]->dip;

						Point2D mh = h * g;
						dip += d + mh;
#ifdef calcOp
						op += 2;
#endif
						if (order >= 2)
						{
							const Point2D& q = mChildren[s]->qua;
							Point2D dh = multz(d, h);
							mh = multz(mh, h);
							qua += q + 2.0 * dh + mh;
#ifdef calcOp
							op += 2;
#endif
							if (order >= 3)
							{
								const Point2D& o = mChildren[s]->oct;

								Point2D qh = multz(q, h);
								dh = multz(dh, h);
								mh = multz(mh, h);
								oct += o + 3.0 * qh + 3.0 * dh + mh;
#ifdef calcOp					
								op += 6;
#endif					
								if (order >= 4)
								{
									const Point2D& hh = mChildren[s]->hex;

									Point2D oh = multz(o, h);
									qh = multz(qh, h);
									dh = multz(dh, h);
									mh = multz(mh, h);

									hex += hh + 4.0 * oh + 6.0 * qh + 4.0 * dh + mh;
#ifdef calcOp						
									op += 8;
#endif
								}//if (order >= 4)
							}//if (order >= 3)
						}//if (order >= 2)
					}
				}
			}
		}//if(!lowLevel)
		else
		{
#ifdef pan
			CalcPanelsParams();
#else
			CalcPointsParams();
#endif
		}//else
	}//CalculateTreeParams()

	// Вычисление параметров нижней ячейки по набору панелей (для решения системы)
	void Tree::CalcPanelsParams()
	{
		mon = 0.0;
		dip.toZero();
		qua.toZero();
		oct.toZero();
		hex.toZero();

		double gDop, gLin;
		Point2D h;
		for (size_t i : index)
		{
			const PointsCopy& itDop = points[i];
			gDop = itDop.g();
			gLin = itDop.gamLin;

			const double g = gDop;
			mon += g;

			h = (itDop.r() - posCentre);

			const Point2D& d = itDop.dipP * gLin; ////
			Point2D mh = h * g;
			dip += d + mh;
#ifdef calcOp
			op += 2;
#endif
			if (order >= 2)
			{
				const Point2D& q = itDop.quaP * gDop;
#ifdef calcOp
				op += 2;
#endif
				Point2D dh = multz(d, h); ////
				mh = multz(mh, h);
				qua += q + dh + 0.5 * mh;
#ifdef calcOp
				op += 2;
#endif
				if (order >= 3)
				{
					Point2D qh = multz(q, h);
					dh = multz(dh, h);		//// 
					mh = multz(mh, h);
					const Point2D& o = itDop.octP * gLin;

					oct += o + 1.5 * qh + 0.75 * dh + 0.25 * mh;
#ifdef calcOp					
					op += 4;
#endif				
					if (order >= 4)
					{
						const Point2D& hh = itDop.hexP * gDop;
#ifdef calcOp
						op += 2;
#endif
						Point2D oh = multz(o, h); //// 
						qh = multz(qh, h);
						dh = multz(dh, h);		//// 
						mh = multz(mh, h);

						hex += hh + 2.0 * oh + 1.5 * qh + 0.5 * dh + 0.125 * mh;
#ifdef calcOp						
						op += 4;
#endif
					}//if (order >= 4)
				}//if (order >= 3)
			}//if (order >= 2)
		}//for i
	}//CalcPanelsParams()

	// Вычисление параметров нижней ячейки по набору вихрей (для расчета скоростей)
	void Tree::CalcPointsParams()
	{
		mon = 0.0;
		dip.toZero();
		qua.toZero();
		oct.toZero();
		hex.toZero();

		double gDop;
		Point2D dr;
		for (size_t i : index)
		{
			const PointsCopy& itDop = points[i];
			gDop = itDop.g();
			mon += gDop;

			dr = (itDop.r() - posCentre);
			Point2D z1 = gDop * dr;
			dip += z1;
#ifdef calcOp
			op += 2;
#endif
			if (order >= 2)
			{

				Point2D z2 = multz(dr, z1);
				qua += z2;
#ifdef calcOp
				op += 2;
#endif
				if (order >= 3)
				{
					Point2D z3 = multz(dr, z2);
					oct += z3;
#ifdef calcOp
					op += 2;
#endif
					if (order >= 4)
					{
						Point2D z4 = multz(dr, z3);
						hex += z4;
#ifdef calcOp
						op += 2;
#endif
					}
				}//order 3
			}//order 2
		}//for i
	}//CalcPointsParams()

	void Tree::ClearCoeff()
	{
		E.toZero({ 0.0, 0.0 });
	}

	//Расчет коэффициентов разложения в ряд Тейлора внутри ячейки нижнего уровня
	void Tree::CalcLocalCoeffToLowLevel(Tree* lowTree, bool calcCloseTrees)
	{
		if (lowTree != this)
		{
			Tree& lt = *lowTree;

			double h, h0;
			h = fabs(this->iRight - this->iLeft) + fabs(this->iTop - this->iBottom);
			h0 = fabs(lt.iRight - lt.iLeft) + fabs(lt.iTop - lt.iBottom);

			Point2D& r0 = lt.posCentre;
			Point2D& r1 = this->posCentre;

			double crit = dist(r0, r1);
#ifdef calcOp
			op += 3;
			op += 2;
#endif
			// если выполнен критерий дальности => считаем коэффициенты
			if ((crit >= (h0 + h + 2.0 * eps) / theta))
			{
				Point2D  rr = r0 - r1;

				double inrm2 = 1.0 / rr.length2();

#ifdef calcOp
				op += 3;
#endif

				Point2D cftr = inrm2 * rr;
				lt.E[0] += mon * cftr;
				Point2D thD = inrm2 * multz(cftr, rr);
				lt.E[0] += multzA(thD, dip);

#ifdef calcOp
				op += 7;
#endif

				if (order >= 2)
				{
					Point2D thQ = (2.0 * inrm2) * multz(thD, rr);
#ifdef calcOp
					op += 3;
#endif
					lt.E[0] += 0.5 * multzA(thQ, qua);

					lt.E[1] -= mon * thD;
#ifdef calcOp
					op += 2;
#endif

					if (order >= 3)
					{
						Point2D thO = (3.0 * inrm2) * multz(thQ, rr);
						lt.E[0] += (1.0 / 6.0) * multzA(thO, oct);
#ifdef calcOp
						op += 5;
#endif
						lt.E[1] -= multzA(thQ, dip);
						lt.E[2] += mon * thQ;
#ifdef calcOp
						op += 2;
#endif

						if (order >= 4)
						{
							Point2D thH = (4.0 * inrm2) * multz(thO, rr);
							lt.E[0] += (1.0 / 24.0) * multzA(thH, hex);
#ifdef calcOp
							op += 5;
#endif						
							lt.E[1] -= 0.5 * multzA(thO, qua);
							lt.E[2] += multzA(thO, dip);
							lt.E[3] -= mon * thO;
#ifdef calcOp
							op += 2;
#endif
						}//if order >= 4
					} //if order >= 3
				}//if order >= 2
			}//if crit
			else // если не выполнен критерий, то рекурсия 
			{
				if (!lowLevel)
				{
					mChildren[0]->CalcLocalCoeffToLowLevel(lowTree, calcCloseTrees);
					mChildren[1]->CalcLocalCoeffToLowLevel(lowTree, calcCloseTrees);
				}
				else if (calcCloseTrees)
					lt.closeTrees.push_back(this);
			}
		}//if (lowTree != this)
		else if (calcCloseTrees)
			closeTrees.push_back(this); //себя тоже добавляем в ближнюю зону 
	}//CalcConvCoeffToLowLevel(...)


	void Tree::CalcConvVeloByBiotSavart()
	{
		for (size_t k = 0; k < closeTrees.size(); ++k)
		{
			//Локальные переменные для цикла
			Point2D velI;
			double dst2eps;

			for (size_t i : index)
			{
				PointsCopy& itDop = points[i];

				velI.toZero();

				const Point2D& posI = itDop.r();

				for (size_t j : closeTrees[k]->index)
				{
					const Point2D& posJ = points[j].r();
					const double& gamJ = points[j].g();
					dst2eps = std::max((posI - posJ).length2(), eps2);
					velI += (gamJ / dst2eps) * (posI - posJ).kcross();
#ifdef calcOp
					op += 3;
#endif
				}//for j

				itDop.veloCopy += velI;
			}//for i 
		}//for k
	}//CalcConvByBiotSavart()




	bool isAfter(const PointsCopy& itI, const PointsCopy& itJ)
	{
		return (itI.panBegin == itJ.panEnd);
	}

	// Вспомогательная функция вычисления угла между векторами
	double Alpha(const Point2D& p, const Point2D& s)
	{
		return atan2(cross3(p, s), p & s);
	}

	// Вспомогательная функция вычисления логарифма отношения норм векторов
	double Lambda(const Point2D& p, const Point2D& s)
	{
		return 0.5 * log((s & s) / (p & p));
	}

	// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
	Point2D Omega(const Point2D& a, const Point2D& b, const Point2D& c)
	{
		return (a & b) * c + (Point2D({ -c[1], c[0] })) * cross3(a, b);
	}

	void Tree::CalcInfluenceFromPanelsToPoints()
	{
		double alpha, lambda;

		//auxillary vectors
		Point2D p, s, di, dj;
		Point2D v00;

		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;

		for (size_t i : index)
		{
			PointsCopy& itDop = points[i];

			velI.toZero();

			const Point2D& posI = itDop.r();
			di = itDop.panEnd - itDop.panBegin;
			double ileni = 1.0 / di.length();
			const Point2D& taui = ileni * di;

			itDop.i00.reserve(int(0.1 * points.size()));

			for (size_t k = 0; k < closeTrees.size(); ++k)
			{
				for (size_t j : closeTrees[k]->index)
				{
					const PointsCopy& itDop2 = points[j];

					const Point2D& posJ = itDop2.r();
					const double& gamJ = itDop2.g();

					if (posI != posJ)
					{
						dj = itDop2.panEnd - itDop2.panBegin;
						double ilenj = 1.0 / dj.length();
						const Point2D& tauj = ilenj * dj;

						p = itDop.r() - itDop2.panEnd;
						s = itDop.r() - itDop2.panBegin;

						alpha = Alpha(p, s);

						lambda = Lambda(p, s);

						v00 = tauj;

						tempVel = ilenj * (alpha * v00 + (lambda * v00).kcross());

						velI += gamJ * tempVel;

						if (itDop.i00.capacity() == itDop.i00.size())
							itDop.i00.reserve(itDop.i00.size() * 2);

						itDop.i00.push_back(tempVel);
					}
				}//for j
			}//for k 
			itDop.veloCopy += velI;
		}//for i			
	}//CalcInfluenceFromPanelsToPoints()

	// Расчет влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
	void Tree::CalcInfluenceFromPanels()
	{
		numvector<double, 3> alpha, lambda;

		//auxillary vectors
		Point2D p1, s1, p2, s2, di, dj, i00, i01, i10, i11;
		numvector<Point2D, 3> v00, v11;
		numvector<Point2D, 2> v01, v10;

		int n = (int)points.size();
		//Локальные переменные для цикла
		Point2D velI, velIlin;

		for (size_t i : index)
		{
			PointsCopy& itDop = points[i];
			velI.toZero();
			velIlin.toZero();

			const Point2D& posI = itDop.r();
			di = itDop.panEnd - itDop.panBegin;
			double dilen = di.length();
			double ileni = 1.0 / dilen;

			const Point2D& taui = ileni * di;

			itDop.i00.reserve(int(0.1 * n));

#ifdef linScheme
			itDop.i01.reserve(int(0.1 * n));
			itDop.i10.reserve(int(0.1 * n));
			itDop.i11.reserve(int(0.1 * n));
#endif
			for (size_t k = 0; k < closeTrees.size(); ++k)
			{
				for (size_t j : closeTrees[k]->index)
				{
					const PointsCopy& itDop2 = points[j];
					const Point2D& posJ = itDop2.r();
					const double& gamJ = itDop2.g();
					const double& gamJlin = itDop2.gamLin;

					if (posI != posJ)
					{
						dj = itDop2.panEnd - itDop2.panBegin;
						double djlen = dj.length();
						double ilenj = 1.0 / djlen;
						const Point2D& tauj = ilenj * dj;

						p1 = itDop.panEnd - itDop2.panEnd;
						s1 = itDop.panEnd - itDop2.panBegin;
						p2 = itDop.panBegin - itDop2.panEnd;
						s2 = itDop.panBegin - itDop2.panBegin;


						alpha = { \
									isAfter(itDop2, itDop) ? 0.0 : Alpha(s2, s1), \
									Alpha(s2, p1), \
									isAfter(itDop, itDop2) ? 0.0 : Alpha(p1, p2) \
						};

						lambda = { \
									isAfter(itDop2, itDop) ? 0.0 : Lambda(s2, s1), \
									Lambda(s2, p1), \
									isAfter(itDop, itDop2) ? 0.0 : Lambda(p1, p2) \
						};

						v00 = { \
								Omega(s1, taui, tauj), \
								- Omega(di, taui, tauj), \
								Omega(p2, taui, tauj) \
						};

						i00 = ilenj * (alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2] \
							+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]).kcross());

#ifdef linScheme
						double s1len2 = s1.length2();

						v01 = {
								0.5 / (djlen) * (((p1 + s1) & tauj) * Omega(s1, taui, tauj) - s1len2 * taui),
								-0.5 * dilen / djlen * Omega(s1 + p2, tauj, tauj)
						};

						i01 = ilenj * ((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]\
							+ (((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * dilen * tauj).kcross());

						v10 = {
								-0.5 / dilen * (((s1 + s2) & taui) * Omega(s1, taui, tauj) - s1len2 * tauj),
								0.5 * djlen / dilen * Omega(s1 + p2, taui, taui)
						};

						i10 = ilenj * ((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1] \
							+ (((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * djlen * taui).kcross());

						v11 = {
							1.0 / (12.0 * dilen * djlen) * \
							(2.0 * (s1 & Omega(s1 - 3.0 * p2, taui, tauj)) * Omega(s1, taui, tauj) - \
									s1len2 * (s1 - 3.0 * p2)) - 0.25 * Omega(s1, taui, tauj),
							-dilen / (12.0 * djlen) * Omega(di, tauj, tauj),
							-djlen / (12.0 * dilen) * Omega(dj, taui, taui)
						};


						i11 = ilenj * ((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]\
							+ ((lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
								+ 1.0 / 12.0 * (djlen * taui + dilen * tauj - 2.0 * Omega(s1, taui, tauj))).kcross());
#endif


						/// \todo не считать снова
						if (isAfter(itDop, itDop2))
						{
							itDop.a = ileni * i00;
#ifdef linScheme
							itDop.a1 = ileni * i11;
#endif
						}
						if (isAfter(itDop2, itDop))
						{
							itDop.c = ileni * i00;
#ifdef linScheme
							itDop.c1 = ileni * i11;
#endif
						}


						if (itDop.i00.capacity() == itDop.i00.size())
							itDop.i00.reserve(itDop.i00.size() * 2);

						itDop.i00.push_back(ileni * i00);

						velI += gamJ * i00;
#ifdef linScheme					
						if (itDop.i01.capacity() == itDop.i01.size())
							itDop.i01.reserve(itDop.i01.size() * 2);

						itDop.i01.push_back(ileni * i01);

						if (itDop.i10.capacity() == itDop.i10.size())
							itDop.i10.reserve(itDop.i10.size() * 2);

						itDop.i10.push_back(ileni * i10);


						if (itDop.i11.capacity() == itDop.i11.size())
							itDop.i11.reserve(itDop.i11.size() * 2);

						itDop.i11.push_back(ileni * i11);
						velI += gamJlin * i01;
						velIlin += gamJ * i10 + gamJlin * i11;
#endif
					}
				}//for j
			}//for k 

			(itDop.veloCopy) += ileni * velI;
#ifdef linScheme
			(itDop.veloCopyLin) += ileni * velIlin;
#endif

		}//for i			
	}//CalcInfluenceFromPanels()

	// Обновление влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
	void Tree::UpdateInfluence()
	{
		Point2D velI, velIlin;

		for (size_t i : index)
		{
			PointsCopy& itDop = points[i];
			velI.toZero();
			velIlin.toZero();

			const Point2D& posI = itDop.r();
			int s = 0;

			for (size_t k = 0; k < closeTrees.size(); ++k)
			{
				for (size_t j : closeTrees[k]->index)
				{
					const PointsCopy& itDop2 = points[j];
					const Point2D& posJ = itDop2.r();
					const double& gamJ = itDop2.g();
					const double& gamJlin = itDop2.gamLin;

					if (posI != posJ)
					{
						velI += gamJ * (itDop.i00[s]);
#ifdef linScheme					
						velI += gamJlin * (itDop.i01[s]);
						velIlin += gamJ * (itDop.i10[s]) + gamJlin * (itDop.i11[s]);
#endif
						s++;
					}
				}//for j
			}//for k
			(itDop.veloCopy) += velI;
#ifdef linScheme	
			(itDop.veloCopyLin) += velIlin;
#endif
		}//for i
	}//UpdateInfluence()

	// Расчет конвективных скоростей внутри дерева нижнего уровня
	void Tree::CalcInfluenceFromVortexFarCells()
	{
		for (size_t i : index)
		{
			PointsCopy& itDop = points[i];
			Point2D deltaPos = itDop.r() - posCentre;
			Point2D v = E[0];
			Point2D vL;
			vL.toZero();

			if (order >= 2)
			{
				const Point2D& tau = itDop.tau;
				const double& len = itDop.len;

				Point2D vLin = multzA(E[1], deltaPos);
				v += vLin;
#ifdef linScheme
				vLin = multzA(E[1], (len / 12.0) * tau);
				vL += vLin;
#endif
				if (order >= 3)
				{
					Point2D r2 = multz(deltaPos, deltaPos);
#ifndef pan
					Point2D vSq = 0.5 * multzA(E[2], r2);
#ifdef calcOp
					op += 2;
#endif
#endif

#ifdef pan
					double len2 = len * len;
					Point2D tau2 = multz(tau, tau);
					Point2D vSq = 0.5 * multzA(E[2], r2 + (len2 / 12.0) * tau2);
#ifdef calcOp
					op += 6;
#endif
#endif
					v += vSq;
#ifdef linScheme
					Point2D rtau = multz(deltaPos, tau);
					vSq = 0.5 * multzA(E[2], (len / 6.0) * rtau);
					vL += vSq;
#endif
					if (order >= 4)
					{
						Point2D r3 = multz(r2, deltaPos);
#ifndef pan
						Point2D vCub = 0.16666666666666667 * multzA(E[3], r3);
#ifdef calcOp
						op += 2;
#endif
#endif
#ifdef pan
						Point2D vCub = 0.16666666666666667 * multzA(E[3], r3 + 0.25 * len2 * multz(tau2, deltaPos));
#ifdef calcOp
						op += 5;
#endif
#endif
						v += vCub;
#ifdef linScheme
						vCub = 0.16666666666666667 * multzA(E[3], 0.25 * len * multz(deltaPos, rtau) + 0.0125 * len2 * len * multz(tau2, tau));
						vL += vCub;
#endif
					}//if order >= 4
				}//if order >= 3
			}//if order >= 2

			itDop.veloCopy += Point2D({ -v[1], v[0] });
			itDop.veloCopy *= IDPI;
#ifdef linScheme
			itDop.veloCopyLin += Point2D({ -vL[1], vL[0] });
			itDop.veloCopyLin *= IDPI;
#endif
		}
	}//CalcConvVeloToLowLevel()


	void Tree::Gauss(std::vector<std::vector<double>>& A, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();

		numvector<double, 3> alpha, lambda;

		//auxillary vectors
		Point2D p1, s1, p2, s2, di, dj, i00, i01, i10, i11;
		numvector<Point2D, 3> v00, v11;
		numvector<Point2D, 2> v01, v10;

		//Локальные переменные для цикла
		Point2D velI, velIlin;

		for (int i = 0; i < n; ++i)
		{
			const PointsCopy& itDop = pnt[i];
			velI.toZero();
			velIlin.toZero();

			const Point2D& posI = itDop.r();
			di = itDop.panEnd - itDop.panBegin;
			double dilen = di.length();
			double ileni = 1.0 / dilen;

			const Point2D& taui = ileni * di;

			for (int j = 0; j < n; ++j)
			{
				const PointsCopy& itDop2 = pnt[j];
				const Point2D& posJ = itDop2.r();

				if (posI != posJ)
				{
					dj = itDop2.panEnd - itDop2.panBegin;
					double djlen = dj.length();
					double ilenj = 1.0 / djlen;
					const Point2D& tauj = ilenj * dj;

					p1 = itDop.panEnd - itDop2.panEnd;
					s1 = itDop.panEnd - itDop2.panBegin;
					p2 = itDop.panBegin - itDop2.panEnd;
					s2 = itDop.panBegin - itDop2.panBegin;


					alpha = { \
								isAfter(itDop2, itDop) ? 0.0 : Alpha(s2, s1), \
								Alpha(s2, p1), \
								isAfter(itDop, itDop2) ? 0.0 : Alpha(p1, p2) \
					};

					lambda = { \
								isAfter(itDop2, itDop) ? 0.0 : Lambda(s2, s1), \
								Lambda(s2, p1), \
								isAfter(itDop, itDop2) ? 0.0 : Lambda(p1, p2) \
					};

					v00 = { \
							Omega(s1, taui, tauj), \
							- Omega(di, taui, tauj), \
							Omega(p2, taui, tauj) \
					};

					i00 = ilenj * (alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2] \
						+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]).kcross());

#ifdef linScheme
					double s1len2 = s1.length2();

					v01 = {
							0.5 / (djlen) * (((p1 + s1) & tauj) * Omega(s1, taui, tauj) - s1len2 * taui),
							-0.5 * dilen / djlen * Omega(s1 + p2, tauj, tauj)
					};

					i01 = ilenj * ((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]\
						+ (((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * dilen * tauj).kcross());

					v10 = {
							-0.5 / dilen * (((s1 + s2) & taui) * Omega(s1, taui, tauj) - s1len2 * tauj),
							0.5 * djlen / dilen * Omega(s1 + p2, taui, taui)
					};

					i10 = ilenj * ((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1] \
						+ (((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * djlen * taui).kcross());

					v11 = {
						1.0 / (12.0 * dilen * djlen) * \
						(2.0 * (s1 & Omega(s1 - 3.0 * p2, taui, tauj)) * Omega(s1, taui, tauj) - \
								s1len2 * (s1 - 3.0 * p2)) - 0.25 * Omega(s1, taui, tauj),
						-dilen / (12.0 * djlen) * Omega(di, tauj, tauj),
						-djlen / (12.0 * dilen) * Omega(dj, taui, taui)
					};


					i11 = ilenj * ((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]\
						+ ((lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
							+ 1.0 / 12.0 * (djlen * taui + dilen * tauj - 2.0 * Omega(s1, taui, tauj))).kcross());
#endif


					A[i][j] = ileni * (i00 & itDop.tau);
					A[i][j + n] = ileni * (i01 & itDop.tau);
					A[i + n][j] = ileni * (i10 & itDop.tau);
					A[i + n][j + n] = ileni * (i11 & itDop.tau);
				}
			}//for j

			A[i][i] = -0.5 * ileni;
			A[i + n][i + n] = -0.041666666666666667 * ileni;

		}
	}

}//namespace BH

