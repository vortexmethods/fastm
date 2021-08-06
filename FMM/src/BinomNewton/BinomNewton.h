/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.0    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2021/08/05     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: BinomNewton.h                                                    |
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
\brief Вычисление биномиальных коэффициентов
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.0
\date 05 августа 2021 г.
*/

#ifndef BINOMNEWTON_H_
#define BINOMNEWTON_H_

#include <vector>

namespace FMM
{

	class BinomNewton
	{
	private:
		std::vector<double> cft;
		int dim;

	public:
		BinomNewton(int n)
			: dim(n + 1)
		{
			cft.resize(dim * dim, 0.0);
			cft[0 * dim + 0] = 1.0;

			for (int i = 1; i <= n; ++i)
			{
				cft[i * dim + 0] = 1.0;
				cft[i * dim + i] = 1.0;
				for (int j = 1; j < i; ++j)
					cft[i * dim + j] = cft[(i - 1) * dim + j] + cft[(i - 1) * dim + (j - 1)];
			}
		}

		double operator()(int p, int q)
		{
			return cft[p * dim + q];
		}

	};

}//namespace FMM

#endif