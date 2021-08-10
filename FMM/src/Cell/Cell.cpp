/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.0    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2021/08/05     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Cell.cpp                                                         |
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
\brief ���������� ������ Cell
\author ���������� ���� ��������������
\author ������ ������� ��������
\author �������� ����� ��������
\version 1.0
\date 05 ������� 2021 �.
*/

#include <algorithm>
#include <iterator>

#include "Cell.h"

namespace FMM
{

	Cell::Cell(std::vector<Particle>& pp_, bool unit) : pp(pp_)
	{
		if (unit)
		{
			r0 = { 0.0, 0.0 };
			wh = { 1.0, 1.0 };
			c = { r0[0] + 0.5 * wh[0], r0[1] + 0.5 * wh[1] };
		}
		else
		{
			auto cmpX = [](const Particle& a, const Particle& b) {return a.r[0] < b.r[0]; };
			auto cmpY = [](const Particle& a, const Particle& b) {return a.r[1] < b.r[1]; };


			auto itMaxX = std::max_element(pp.begin(), pp.end(), cmpX);
			auto itMaxY = std::max_element(pp.begin(), pp.end(), cmpY);

			auto itMinX = std::min_element(pp.begin(), pp.end(), cmpX);
			auto itMinY = std::min_element(pp.begin(), pp.end(), cmpY);

			r0 = { itMinX->r[0], itMinY->r[1] };
			wh = { itMaxX->r[0], itMaxY->r[1] };

			wh -= r0;

			r0 -= 0.01 * wh;
			wh *= 1.02;
			c = { r0[0] + 0.5 * wh[0], r0[1] + 0.5 * wh[1] };

			double maxDim = std::max(wh[0], wh[1]);
			wh = { maxDim, maxDim };
			r0 = { c[0] - 0.5 * maxDim, c[1] - 0.5 * maxDim };
		}


		particles.reserve(pp.size());
		for (int i = 0; i < (int)pp.size(); ++i)
			//particles.push_back(pp.begin() + i);
			particles.push_back(i);

		leaf = false;
		level = 0;
		parent = -1;
		index = 0;

	}


	Cell::Cell(const Cell& parent_, int index_) //�������
		: pp(parent_.pp)
	{
		const Point2D& pr0 = parent_.r0;
		const Point2D& pwh = parent_.wh;

		switch (index_)
		{
		case 0:
			r0 = { pr0[0], pr0[1] + 0.5 * pwh[1] };
			wh = { pr0[0] + 0.5 * pwh[0], pr0[1] + pwh[1] };
			break;

		case 1:
			r0 = { pr0[0] + 0.5 * pwh[0], pr0[1] + 0.5 * pwh[1] };
			wh = { pr0[0] + pwh[0], pr0[1] + pwh[1] };
			break;

		case 2:
			r0 = { pr0[0] , pr0[1] };
			wh = { pr0[0] + 0.5 * pwh[0], pr0[1] + 0.5 * pwh[1] };
			break;

		case 3:
			r0 = { pr0[0] + 0.5 * pwh[0], pr0[1] };
			wh = { pr0[0] + pwh[0], pr0[1] + 0.5 * pwh[1] };
			break;
		}

		wh -= r0;

		c = { r0[0] + 0.5 * wh[0], r0[1] + 0.5 * wh[1] };

		particles.reserve(parent_.particles.size());

		std::copy_if(parent_.particles.begin(), parent_.particles.end(), back_inserter(particles),
			[&](int i) { return ((pp[i].r[0] <= r0[0] + wh[0]) && (pp[i].r[0] > r0[0]) && (pp[i].r[1] <= r0[1] + wh[1]) && (pp[i].r[1] > r0[1]));  });

		particles.shrink_to_fit();
		level = parent_.level + 1;

		leaf = ((level >= maxLevel) || (particles.size() <= 1)) ? true : false;

		index = index_;
	}


	Cell::~Cell()
	{
	}

}//namespave FMM
