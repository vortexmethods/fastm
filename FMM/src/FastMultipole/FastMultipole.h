/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: FastMultipole.h                                                  |
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
\brief Описание класса FastMultipole
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/


#ifndef FASTMULTIPOLE_H_
#define FASTMULTIPOLE_H_

#include <complex>
#include <fstream>
#include <vector>

#include "BinomNewton.h"
#include "Cell.h"
#include "Particle.h"
#include "Point2D.h"
#include "Tree.h"

namespace FMM
{

	struct FastMultipole
	{
	public:
		static double tDn, tUp, tLeaf, tL2L, tM2L;
		long long opMult, opDiv;
		BinomNewton binom;

		FastMultipole(int n);

		/// Вычисление мультипольных моментов ячеек дерева при обходе от листьев к корню
		// для нижнего уровня считаем по набору вихрей, для остальных --- по потомкам
		void Upward(Tree& qtree, int pos);


		void Downward(Tree& qtree, int pos);



		void M2M(const std::complex<double>* coeffs, Point2D c, std::complex<double>* shift);

		void L2L(const std::complex<double>* coeffs, Point2D c, std::complex<double>* shift);

		void M2L(const std::complex<double>* coeffs, Point2D c, std::complex<double>* inner);

		void potentialDS(Cell& n1, const Cell& n2);

		void forceDS(Cell& n1, const Cell& n2);

		void DirectCalc(Tree& qtree);

		void multipole(Cell& node);

		void InfluenceComputation(Tree& qtree);

		~FastMultipole();

		template <typename T>
		T myPow(T base, unsigned int exp)
		{
			T res = 1;
			while (exp) {
				if (exp & 1)
				{
					res *= base;
#ifdef calcOp
					opMult++;
#endif
				}
				exp >>= 1;
				base *= base;
#ifdef calcOp
				opMult++;
#endif
			}
			return res;
		}
	};


	inline double m1(int n)
	{
		return (n % 2) ? -1.0 : 1.0;
	}


	inline void divComp(const std::complex<double>& a, const std::complex<double>& b, std::complex<double>& c)
	{
		const double& ax = a.real();
		const double& bx = b.real();
		const double& ay = a.imag();
		const double& by = b.imag();

		double zn = bx * bx + by * by;
		c = { (ax * bx + ay * by) / zn, (ay * bx - ax * by) / zn };
	}

	inline void divComp(double& a, const std::complex<double>& b, std::complex<double>& c)
	{
		const double& bx = b.real();
		const double& by = b.imag();

		double zn = bx * bx + by * by;
		c = { (a * bx) / zn, -(a * by) / zn };
	}

}//namespace FMM

#endif