/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: FastMultipole.cpp                                                |
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
\brief Реализация класса FastMultipole
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/


#include <algorithm>
#include <numeric>

#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include <omp.h>

#include "BinomNewton.h"
#include "Cell.h"
#include "FastMultipole.h"
#include "Particle.h"
#include "Point2D.h"
#include "Tree.h"

double mt = 0.0, m2mt = 0.0;

using namespace std::literals;

namespace FMM
{

	double FastMultipole::tDn = 0.0, FastMultipole::tUp = 0.0, FastMultipole::tLeaf = 0.0, FastMultipole::tL2L = 0.0, FastMultipole::tM2L = 0.0;

	FastMultipole::FastMultipole(int n)
		: binom(2 * n)
	{
#ifdef calcOp
		opMult = 0; opDiv = 0;
#endif

	}

	void FastMultipole::InfluenceComputation(Tree& qtree)
	{
		double t1, t2;

		t1 = omp_get_wtime();
		Upward(qtree, 0);
		t2 = omp_get_wtime();

		tUp = t2 - t1;

		for (int i = 0; i < nt; ++i)
			qtree.node[0].inner[i] = 0.0;

		t1 = omp_get_wtime();
		Downward(qtree, 0);

		DirectCalc(qtree);

		t2 = omp_get_wtime();

		tDn = t2 - t1;
#ifdef calcOp
		//std::cout << "opMult = " << 1.0 * (opMult + opDiv + 2 * n) << std::endl;
#endif
	}

	void FastMultipole::DirectCalc(Tree& qtree)
	{
		double t1, t2;
		t1 = omp_get_wtime();

		std::vector<std::vector<Cell>::iterator> leafs;

		for (auto it = qtree.node.begin(); it != qtree.node.end(); ++it)
			if (it->leaf)
				leafs.push_back(it);

		std::complex<double> diff, z;

#ifndef calcOp
#pragma omp parallel for schedule(guided,4) private(diff,z)
#endif
		for (int i = 0; i < (int)leafs.size(); ++i)
		{
			Cell& leaf = *leafs[i];

			std::complex<double> pows[nt];
			pows[0] = 1.0;

			for (int j = 0; j < (int)leaf.particles.size(); j++)
			{
				diff = leaf.pp[leaf.particles[j]].z - std::complex<double>{leaf.c[0], leaf.c[1]};

				//Если потенциал не нужен - суммируем до nt-1, иначе - до nt
				for (int q = 1; q < (calcPot ? nt : nt - 1); ++q)
				{
					pows[q] = pows[q - 1] * diff;
#ifdef calcOp 
					opMult += 4;
#endif
				}

				if (calcPot)
					z = 0.0;

				//Если потенциал не нужен - суммируем с единицы, иначе - с нуля
				for (int k = (calcPot ? 0 : 1); k < nt; k++)
				{
					if (calcPot)
					{
						z += leaf.inner[k] * pows[k];
#ifdef calcOp
						opMult += 4;
#endif
					}
					leaf.pp[leaf.particles[j]].f += (double)k * leaf.inner[k] * pows[k - 1];
#ifdef calcOp
					opMult += 6;
#endif
				}

				if (calcPot)
					leaf.pp[leaf.particles[j]].phi = real(z);
			}


			for (int j = 0; j < (int)leaf.closeNeighbor.size(); j++)
			{
				//if (calcPot)
					potentialDS(leaf, qtree.node[leaf.closeNeighbor[j]]);

				forceDS(leaf, qtree.node[leaf.closeNeighbor[j]]);
			}

			//if (calcPot)
				potentialDS(leaf, leaf);

			forceDS(leaf, leaf);
		}

		t2 = omp_get_wtime();
		tLeaf += t2 - t1;
	}

	void FastMultipole::Upward(Tree& qtree, int pos)
	{
		for (int i = 0; i < nt; ++i)
			qtree.node[pos].inner[i] = 0.0;

		if ((qtree.node[pos].leaf))
		{
			double t1 = omp_get_wtime();
			multipole(qtree.node[pos]);
			mt += omp_get_wtime() - t1;
		}
		else
		{
			std::complex<double> temp[nt];
			Point2D c = {};
			for (int i = 0; i < (int)qtree.node[pos].child.size(); i++)
			{
				Upward(qtree, qtree.node[pos].child[i]);
				c = qtree.node[qtree.node[pos].child[i]].c - qtree.node[pos].c;
				double t1 = omp_get_wtime();
				M2M(qtree.node[qtree.node[pos].child[i]].outer, c, temp);
				m2mt += omp_get_wtime() - t1;
				for (int j = 0; j < nt; j++)
				{
					qtree.node[pos].outer[j] += temp[j];
				}
			}
		}
	}

	void FastMultipole::Downward(Tree& qtree, int pos)
	{
		Point2D c;
		std::complex<double> temp[nt], diff, z;

		double t1, t2;

		for (int i = 0; i < (int)qtree.node[pos].farNeighbor.size(); ++i)
		{
			const auto& nd = qtree.node[qtree.node[pos].farNeighbor[i]];

			c = nd.c - qtree.node[pos].c;

			t1 = omp_get_wtime();

			M2L(nd.outer, c, temp);

			t2 = omp_get_wtime();
			tM2L += t2 - t1;

			for (int j = 0; j < nt; j++)
				qtree.node[pos].inner[j] += temp[j];

		}

		if (!qtree.node[pos].leaf)
		{
			for (int i = 0; i < (int)qtree.node[pos].child.size(); i++)
			{
				const int& ch = qtree.node[pos].child[i];

				c = qtree.node[pos].c - qtree.node[ch].c;


				t1 = omp_get_wtime();
				L2L(qtree.node[pos].inner, c, qtree.node[ch].inner);
				t2 = omp_get_wtime();
				tL2L += t2 - t1;

				Downward(qtree, ch);
			}
		}
	}

	void FastMultipole::L2L(const std::complex<double>* coeffs, Point2D c, std::complex<double>* shift)
	{
		std::complex <double> z0 = { c[0],c[1] };
/*
		for (int i = 0; i < nt; ++i)
		{
			shift[i] = 0.0;
			for (int j = i; j < nt; ++j)
			{
				shift[i] += coeffs[j] * binom(j, i) * myPow(-z0, (j - i));
#ifdef calcOp
				opMult += 6;
#endif
			}
		}
*/
		for (int j = 0; j < nt; ++j)
			shift[j] = coeffs[j];

		for (int j = 0; j < nt; ++j)
		{			
			for (int k = nt - j - 1; k < nt - 1; ++k)
			{
				shift[k] -= z0 * shift[k + 1];
			}
		}


		//return shift;
	}

	double logCompRe(const std::complex<double>& a)
	{
		return 0.5 * log(a.real() * a.real() + a.imag() * a.imag());
	}

	std::complex<double> iz0powi[nt], prod[nt];

	void FastMultipole::M2L(const std::complex<double>* coeffs, Point2D c, std::complex<double>* inner)
	{		
		std::complex<double> z0 = { c[0], c[1] };

		iz0powi[0] = 1.0;
		prod[0] = 0.0;

		std::complex<double> tmp;

		if (calcPot)
		{
			inner[0] = coeffs[0] * 0.5 * log(c[0] * c[0] + c[1] * c[1]);
#ifdef calcOp
			opMult += 6;
#endif
		}

		for (int i = 1; i < nt; ++i)
		{
			//iz0powi[i] = iz0powi[i - 1] / z0;
			divComp(iz0powi[i - 1], z0, iz0powi[i]);


			prod[i] = m1(i) * coeffs[i] * iz0powi[i];
#ifdef calcOp
			opMult += 2; opMult += 4;
			opMult += 6; opDiv += 2;
#endif
			inner[0] += prod[i];
		}

		for (int i = 1; i < nt; ++i)
		{
			inner[i] = (-1.0 / i) * coeffs[0] * iz0powi[i];
#ifdef calcOp
			opDiv += 1; opMult += 6;
#endif
			tmp = 0.0;
			for (int j = 1; j < nt; ++j)
			{
				tmp += binom(i + j - 1, j - 1) * prod[j];
#ifdef calcOp
				opMult += 2;
#endif
			}

			inner[i] += iz0powi[i] * tmp;
#ifdef calcOp
			opMult += 2;
#endif

		}
	}

	void FastMultipole::forceDS(Cell& n1, const Cell& n2)
	{
		double xi, yi, xj, yj;
		double znam;

		std::complex<double> res;
		double r2 = 0.0;
		for (int i = 0; i < (int)n1.particles.size(); ++i)
		{
			xi = reinterpret_cast<double(&)[2]>(n1.pp[n1.particles[i]].z)[0];
			yi = reinterpret_cast<double(&)[2]>(n1.pp[n1.particles[i]].z)[1];

			for (int j = 0; j < (int)n2.particles.size(); ++j)
			{
				xj = reinterpret_cast<double(&)[2]>(n2.pp[n2.particles[j]].z)[0];
				yj = reinterpret_cast<double(&)[2]>(n2.pp[n2.particles[j]].z)[1];

				znam = n2.pp[n2.particles[j]].q / std::max((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj), eps2);

				n1.pp[n1.particles[i]].f += std::complex<double>{(xi - xj)* znam, -(yi - yj) * znam};
#ifdef calcOp
				opMult += 2;
				opMult += 2; opDiv += 1;
#endif
			}
		}
	}

	void FastMultipole::potentialDS(Cell& n1, const Cell& n2)
	{
		double r2 = 0.0;
		for (int i = 0; i < (int)n1.particles.size(); i++)
		{
			for (int j = 0; j < (int)n2.particles.size(); j++)
			{
				r2 = n1.pp[n1.particles[i]].r.dist2To(n2.pp[n2.particles[j]].r);
#ifdef calcOp
				opMult += 2;
#endif
				if (r2 > 1e-12)
				{
					n1.pp[n1.particles[i]].phi += n2.pp[n2.particles[j]].q * 0.5 * log(r2);
#ifdef calcOp
					opMult += 2;
#endif
				}
			}
		}
	}

	std::vector<std::complex<double>> tabpow(nt, 1.0);
	
	void FastMultipole::M2M(const std::complex<double>* coeffs, Point2D c, std::complex<double>* shift)
	{
		std::complex <double> z0 = { c[0], c[1] };
		
		shift[0] = coeffs[0];

		
		for (int i = 1; i < nt; ++i)
			tabpow[i] = tabpow[i - 1] * z0;

		for (int i = 1; i < nt; ++i)
		{
			shift[i] = -(coeffs[0] * tabpow[i]) / (double)i;
#ifdef calcOp
			opMult += 4; opDiv += 2;
#endif		
			for (int j = 1; j <= i; ++j)
			{
				//double b = binom(i - 1, j - 1);
				shift[i] += coeffs[j] * tabpow[i-j] * binom(i - 1, j - 1);
#ifdef calcOp
				opMult += 4; opMult += 2;
#endif
			}
		}
	}



	void FastMultipole::multipole(Cell& node)
	{
/*		
		std::complex <double> z0 = {node.c[0], node.c[1]};
		//int nd;
		for (int i = 0; i < (int)node.particles.size(); i++)
		{
			const auto& nd = node.pp[node.particles[i]];
			//nd = node.particles[i];

			node.outer[0] += nd.q; //node.pp[nd].q;

			for (int j = 1; j < nt; j++)
			{
				//node.outer[j] -= node.pp[nd].q / j * myPow(node.pp[nd].z - z0, j);
				node.outer[j] -= nd.q / j * myPow(nd.z - z0, j);
#ifdef calcOp
					opDiv += 1; opMult += 2;
#endif
			}
		}
//*/

		//*		
		std::complex <double> z0 = {node.c[0], node.c[1]};
		std::complex <double> tmp;

		for (int j = 1; j < nt; ++j)
		{
			//tmp = 0.0;
			for (int i = 0; i < (int)node.particles.size(); ++i)
			{
				const auto& nd = node.pp[node.particles[i]];
				node.outer[j] -= nd.q * myPow(nd.z - z0, j);
			}			
			node.outer[j] /= (double)j;
		}//for j

		for (int i = 0; i < (int)node.particles.size(); i++)
		{
			const auto& nd = node.pp[node.particles[i]];
			node.outer[0] += nd.q;
		}
		//*/

		/*
		std::complex <double> z0 = { node.c[0], node.c[1] };
		for (int i = 1; i < nt; ++i)
		{
			auto temp = std::accumulate(node.particles.begin(), node.particles.end(), 0.0i, [&z0, i, &node, this](const std::complex<double>& x, const int& y)
				{
					return x + node.pp[y].q * myPow(node.pp[y].z - z0, i);
				});
			node.outer[i] = -temp / double(i);
		}

		node.outer[0] = accumulate(node.particles.begin(), node.particles.end(), 0.0i, [this, &node](const std::complex<double>& x, const int& y) { return x + node.pp[y].q; });
		*/

	}

	FastMultipole::~FastMultipole()
	{
	}

}//namespace FMM