/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Tree.cpp                                                         |
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
\brief Реализация класса Tree
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include <omp.h>

#include "Tree.h"

namespace FMM
{

	Tree::Tree(std::vector<Particle>& pp)
	{
		double ta = omp_get_wtime();

		node.emplace_back(Cell(pp, false));

		buildTree(node.back(), (int)node.size() - 1);

		double tb = omp_get_wtime();

		//for (int c = 0; c < (int)node.size(); ++c)
		//{
		//	node[c].closeNeighbor.reserve(8);
		//	fillNeighbors(node[c].closeNeighbor, c, 0);
		//}

		for (int c = 1; c < (int)node.size(); ++c)
		{
			node[c].closeNeighbor.reserve(8);
			Cell p = node[node[c].parent];
			for (int j = 0; j < (int)p.child.size(); ++j)
			{
				if (p.child[j] != c)
				{
					node[c].closeNeighbor.push_back(p.child[j]);
					//node[p.child[j]].closeNeighbor.push_back(c);
				}
			}
			for (int i = 0; i < (int)p.closeNeighbor.size(); ++i)
			{
				fillNeighbors(node[c].closeNeighbor, c, p.closeNeighbor[i]);
				Cell cN = node[p.closeNeighbor[i]];
				for (int j = 0; j < (int)cN.child.size(); ++j)
					fillNeighbors(node[c].closeNeighbor, c, cN.child[j]);
			}

			//fillNeighbors(node[c].closeNeighbor, c, 0);
		}

		double tc = omp_get_wtime();

		for (int c = 1; c < (int)node.size(); ++c)
			if (node[c].level > 1) farNeighbors(c);

		double td = omp_get_wtime();

		//std::cout << "build = " << tb - ta << ", neib. = " << tc - tb << ", neibfar. = " << td - tc << std::endl;

		//std::cout << "full = " << td - ta << endl;

		//fileTree();
	}

	void Tree::buildTree(Cell root, int posRoot)
	{
		if (root.leaf)
			return;

		int sz = (int)node.size(), q = 0;
		for (int i = 0; i < 4; ++i)
		{
			//node.emplace_back(Cell(root, i));
			Cell C(root, i);


			node.push_back(C);

			node.back().parent = posRoot;
			++q;
			//node[posRoot].child.emplace_back(node.size() - 1);
			node[posRoot].child.push_back((int)node.size() - 1);

		}

		for (int i = 0; i < 4; ++i)
			buildTree(*(node.begin() + sz + i), sz + i);
	}

	void Tree::fillNeighbors(std::vector<int>& found, int who, int trav)
	{
		double EPS = 1e-10;

		if ((node[trav].leaf) || (node[trav].level == node[who].level))
		{
			if ((fabs(fabs(node[who].c[0] - node[trav].c[0]) - EPS) < 0.5 * (node[who].wh[0] + node[trav].wh[0])) && (fabs(fabs(node[who].c[1] - node[trav].c[1]) - EPS) < 0.5 * (node[who].wh[1] + node[trav].wh[1])))
			{
				found.push_back(trav);
			}
			return;
		}
		return;
	}

	void Tree::farNeighbors(int who)
	{
		node[who].farNeighbor.reserve(32);

		std::set<int> neibst(node[who].closeNeighbor.begin(), node[who].closeNeighbor.end());
		neibst.insert(who);

		std::set<int> st;

		//Цикл по ближайшим соседям
		for (int i = 0; i < (int)node[who].closeNeighbor.size(); ++i)
		{
			const Cell& nb = node[node[who].closeNeighbor[i]];

			//если сосед того же уровня, что ячейка
			if (nb.level == node[who].level)
			{
				//Цикл по соседям соседа
				for (int j = 0; j < (int)nb.closeNeighbor.size(); ++j)
				{
					const Cell& nbnb = node[nb.closeNeighbor[j]];

					if (nbnb.level == node[who].level)
					{
						for (int k = 0; k < (int)node[nbnb.parent].child.size(); ++k)
							st.insert(node[nbnb.parent].child[k]);
					}
					else st.insert(nb.closeNeighbor[j]);
				}
			}
		}
		std::set_difference(st.begin(), st.end(), neibst.begin(), neibst.end(), std::back_inserter(node[who].farNeighbor));
		node[who].farNeighbor.shrink_to_fit();
	}

	void Tree::fileTree()
	{
		std::ofstream fin1("tree.txt");
		std::ofstream fin2("particles.txt");
		std::ofstream fin3("q.txt");
		for (int i = 0; i < (int)node.size(); i++)
		{
			//fin << /*"{" <<*/qtree[i].r0 << " , " << qtree[i].r0 + qtree[i].wh /*<< "}"*/ <<endl;
			fin1 << node[i].r0 << std::endl;
			fin1 << node[i].wh + node[i].r0 << std::endl;

		}
		for (int i = 0; i < (int)node[0].pp.size(); i++)
		{
			fin2 << node[0].pp[i] << std::endl;
			fin3 << node[0].pp[i].q << std::endl;
		}
		fin1.close();
		fin2.close();
		fin3.close();
	}

	Tree::~Tree()
	{
	}

}//namespace FMM
