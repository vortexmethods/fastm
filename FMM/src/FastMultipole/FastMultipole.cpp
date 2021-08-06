/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.0    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2021/08/05     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
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
\version 1.0
\date 05 августа 2021 г.
*/


#include <algorithm>
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
				if (calcPot)
					potentialDS(leaf, qtree.node[leaf.closeNeighbor[j]]);

				forceDS(leaf, qtree.node[leaf.closeNeighbor[j]]);
			}

			if (calcPot)
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
			multipole(qtree.node[pos]);
		else
		{
			std::complex<double> temp[nt];
			Point2D c = {};
			for (int i = 0; i < (int)qtree.node[pos].child.size(); i++)
			{
				Upward(qtree, qtree.node[pos].child[i]);
				c = qtree.node[qtree.node[pos].child[i]].c - qtree.node[pos].c;
				M2M(qtree.node[qtree.node[pos].child[i]].outer, c, temp);
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
		//return shift;
	}

	double logCompRe(const std::complex<double>& a)
	{
		return 0.5 * log(a.real() * a.real() + a.imag() * a.imag());
	}

	void FastMultipole::M2L(const std::complex<double>* coeffs, Point2D c, std::complex<double>* inner)
	{
		std::complex<double> iz0powi[nt], prod[nt];
		std::complex<double> z0 = { c[0], c[1] };

		iz0powi[0] = 1.0;
		prod[0] = 0.0;

		if (calcPot)
		{
			inner[0] = coeffs[0] * 0.5 * log(c[0] * c[0] + c[1] * c[1]);
#ifdef calcOp
			opMult += 6;
#endif
		}

		for (int i = 1; i < nt; ++i)
		{
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
			for (int j = 1; j < nt; ++j)
			{
				inner[i] += binom(i + j - 1, j - 1) * prod[j] * iz0powi[i];
#ifdef calcOp
				opMult += 6;
#endif
			}
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

	void FastMultipole::M2M(const std::complex<double>* coeffs, Point2D c, std::complex<double>* shift)
	{
		std::complex <double> z0 = { c[0], c[1] };
		for (int i = 0; i < nt; i++)
		{
			if (i == 0)
			{
				shift[i] = coeffs[0];
				continue;
			}

			shift[i] = -(coeffs[0] * myPow(z0, i)) / (double)i;
#ifdef calcOp
			opMult += 4; opDiv += 2;
#endif		
			for (int j = 1; j <= i; j++)
			{
				//double b = binom(i - 1, j - 1);
				shift[i] += coeffs[j] * myPow(z0, (i - j)) * binom(i - 1, j - 1);
#ifdef calcOp
				opMult += 4; opMult += 2;
#endif
			}
		}
	}

	void FastMultipole::multipole(Cell& node)
	{
		std::complex <double> z0 = { node.c[0], node.c[1] };
		int nd;
		for (int i = 0; i < (int)node.particles.size(); i++)
		{
			nd = node.particles[i];
			for (int j = 0; j < nt; j++)
			{
				node.outer[j] += (j == 0) ? node.pp[nd].q : (-node.pp[nd].q / j * myPow(node.pp[nd].z - z0, j));

				if (j != 0)
				{
#ifdef calcOp
					opDiv += 1; opMult += 2;
#endif
				}
			}
		}

	}

	FastMultipole::~FastMultipole()
	{
	}

}//namespace FMM