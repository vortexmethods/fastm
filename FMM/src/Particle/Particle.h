/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.0    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2021/08/05     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Particle.h                                                       |
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
\brief Описание структуры Particle
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.0
\date 05 августа 2021 г.
*/



#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <cmath>
#include <complex>
#include <ostream>

#include "Point2D.h"

namespace FMM
{

	struct Particle
	{
	public:
		std::complex <double> z, f;
		Point2D r;
		double q, phi;

		bool operator==(const Particle& p2) const
		{
			return (r == p2.r);
		}

		bool operator!=(const Particle& p2) const
		{
			return (r != p2.r);
		}

		bool operator<(const Particle& p2) const
		{
			return (r[0] < p2.r[0]) || ((r[0] == p2.r[0]) && (r[1] < p2.r[1]));
		}

		static int counter;

		static Particle getParticle()
		{
			/*Point2D r = { 0.5 + 0.499*cos(2.0*3.1415926535*(counter + 0.5) / 50.0), 0.5 + 0.499*sin(2.0*3.1415926535*(counter + 0.5) / 50.0) };
			double phi = 0.0;
			complex<double> f = 0.0;
			double q = 1.0;
			complex <double> z = { r[0], r[1] };
			++counter;

			return { z, f, r, q, phi };*/


			Point2D r = { (rand() % 10000) * 0.0001 + 0.00004, (rand() % 10000) * 0.0001 + 0.00004 };
			double phi = 0.0;
			std::complex<double> f = 0.0;
			double q = (rand() % 100) * 0.01;
			std::complex <double> z = { r[0], r[1] };
			return { z, f, r, q, phi };
		};

		~Particle() {};
	};

	std::ostream& operator << (std::ostream& str, const Particle& p);

}//namespace FMM

#endif