/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: PointsCopy.h                                                     |
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
\brief Заголовок вспомогательного класса
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/


#ifndef POINTSCOPY_H_
#define POINTSCOPY_H_

#include "Vortex2D.h"
#include "Params.h"
#include "omp.h"

namespace BH
{

	class PointsCopy : public Vortex2D
	{
	public:
		Point2D veloCopy, veloCopyLin;
		double velTau, velTauLin;
		std::vector<Point2D> i00;
		Point2D a, c; //для предобуславливателя в Т0

#ifdef linScheme
		std::vector<Point2D> i01, i10, i11;
		Point2D a1, c1; //для предобуславливателя в Т1
		double gamLin;
#endif

		Point2D panBegin, panEnd;
		std::vector<Point2D> Mpan;

#ifdef asympScheme
		std::vector<Point2D> MpanAs;
#endif



		Point2D tau;
		double len;

		PointsCopy(const Vortex2D& vtx_) :
			Vortex2D(vtx_), veloCopy({ 0.0, 0.0 }), veloCopyLin({ 0.0, 0.0 }), 
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 })
		{};

		PointsCopy() :
			Vortex2D({ 0.0, 0.0 }, 0.0),
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 }) 
		{};

		PointsCopy(const Point2D& r) :
			Vortex2D(r, 0.0), 
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 }) 
		{};

		~PointsCopy() {};
	};

}//namespace BH

#endif
