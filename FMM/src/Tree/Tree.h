/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.0    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2021/08/05     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Tree.h                                                           |
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
\brief �������� ��������� Tree
\author ���������� ���� ��������������
\author ������ ������� ��������
\author �������� ����� ��������
\version 1.0
\date 05 ������� 2021 �.
*/


#ifndef TREE_H_
#define TREE_H_

#include <vector>

#include "Particle.h"
#include "Cell.h"
#include "Point2D.h"

namespace FMM
{

	struct Tree
	{
	public:

		std::vector<Cell> node;

		Tree(std::vector<Particle>& pp);

		void buildTree(Cell root, int posRoot);

		void fillNeighbors(std::vector<int>& found, int who, int trav);

		void farNeighbors(int who);

		void fileTree();

		~Tree();
	};

}//namespace FMM

#endif