/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Cell.h                                                           |
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
\brief Описание класса Cell
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/


#ifndef CELL_H_
#define CELL_H_

#include <vector>

#include "Params.h"
#include "Particle.h"
#include "Point2D.h"


namespace FMM
{

	struct Cell
	{
	public:
		std::vector<Particle>& pp;
		std::complex <double> outer[nt];
		std::complex <double> inner[nt];

		Point2D wh;
		Point2D r0;
		Point2D c;

		std::vector<int> particles;
		std::vector<int> child;

		bool leaf;
		int level;
		int index;
		int parent;

		std::vector<int> closeNeighbor;
		std::vector<int> farNeighbor;

		Cell(std::vector<Particle>& pp, bool unit = false); //корень
		Cell(const Cell& parent, int index); //потомок
		~Cell();
	};

}//namespace FMM

#endif