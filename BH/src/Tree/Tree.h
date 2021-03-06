/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.0    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2021/08/05     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: Tree.h                                                           |
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
\brief ?????????????????? ???????????? Tree
\author ???????????????????? ???????? ????????????????????????????
\author ???????????? ?????????????? ????????????????
\version 1.0
\date 05 ?????????????? 2021 ??.
*/

#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <memory>

#include "Params.h"
#include "PointsCopy.h"

namespace BH
{

	extern long long op;

	class Tree
	{
	private:
		/// ???????????????? ????????????????????????????
		double iLeft, iRight, iTop, iBottom;

		/// ?????????? ????????????????????????????
		Point2D posCentre;

		/// ?????????? ?????????? ???????????? (?????????????? ????????????)
		int level;
		/// ?????????? ?????????????? ???????????? ????????????
		int p;

		/// ?????????????? ???? ?? ???????????? ???????????????
		bool lowLevel = false;

		/// ?????????????? ?????????????? ???????????? ????????????
		std::unique_ptr<Tree> mChildren[2];

		/// ?????????????????????? ???????????? ???????????? 
		double mon;
		/// ??????????????????, ??????????????????????????, ?????????????????????? ?? ???????????????????????????????? ?????????????? ???????????? 
		Point2D dip, qua, oct, hex;

		/// \brief ???????????????????????? ?????? ???????????????????? ??????????????????
		numvector<Point2D, 4> E;

	public:

		std::vector<PointsCopy>& points;

		/// ???????????? ???????????????????? ???? ?????????????? ?? ?????????????? ???????? (??????, ?????? ???????? ?????????????? ?????????????? "????????????????") 
		//?????????? ?????????? ???????????? ?????? ???????????? ??????????????
		std::vector<Tree*> closeTrees;

		/// ???????????? ???????????????????? ???? ???????????? ?????????????? ???????????? 
		//?????????? ?????????? ???????????? ?????? ??????????
		std::vector<Tree*> lowTrees;

		/// ?????????????? ?????????????????? ?? ???????????????? ?????????????? ??????????
		std::vector<int> index;

		/// ?????????????????? ???? ???????????? (?????? ???????????????????? lowTrees)
		Tree* BigTree;

		std::vector<std::unique_ptr<Tree>> treePtr;

		/// ??????????????????????
		/// \param[in]

		Tree(std::vector<PointsCopy>& points_)
			:points(points_)
		{
			mon = 0.0;
			dip.toZero();
			qua.toZero();
			oct.toZero();
			hex.toZero();

			posCentre.toZero();

			mChildren[0] = nullptr;
			mChildren[1] = nullptr;

			level = 0;

			E.toZero({ 0.0, 0.0 });
		}

		/// ????????????????????
		~Tree() {};


		/// ???????????????????? ?????????? ???????????? ???? ???????????? ???????????????? ????????????
		///
		/// \param[in] points ???????????? ???? ???????????? ????????????, ???? ???????????? ?????????????? ???????????????? ??????????????????????????
		void makeTree(std::vector<PointsCopy>& points);

		/// ???????????????????? ?????????????????? ????????????
		void CreateTree();

		/// ???????????????????? ???????????????????? ???????????? (???????????????????? ?? ?????????????? ??????????????????????????)	
		// ?????? ?????????????? ???????????? ?????????????? ???? ???????????? ????????????, ?????? ?????????????????? --- ???? ????????????????
		void CalculateTreeParams();
		/// ???????????????????? ???????????????????? ???????????? ???????????? ???? ???????????? ???????????? (?????? ?????????????? ??????????????????)
		void CalcPointsParams();
		/// ???????????????????? ???????????????????? ???????????? ???????????? ???? ???????????? ?????????????? (?????? ?????????????? ??????????????)
		void CalcPanelsParams();

		/// ???????????? ???????????????????????? ?????????????????? ???????????? ???????????? ???????????? ???????????? ?????????????? ???????????? ???? ???????????? ???????? ?????????????? ??????????????
		void CalcConvVeloByBiotSavart();

		/// ?????????????????? ?????????????????????????? ?????????????????? ????????????????????
		void ClearCoeff();

		/// ???????????? ?????????????????????????? ???????????????????? ?? ?????? ?????????????? ???????????? ???????????? ?????????????? ????????????
		///
		/// \param[in] lowTree ?????????????????? ???? ???????????? ?????????????? ??????????????
		/// \param[in] calcCloseTrees = true (???? ??????????????????), ???????? ?????????? ?????????????????? ???????????? ?????????????? ??????????
		void CalcLocalCoeffToLowLevel(Tree* lowTree, bool calcCloseTrees = true);

		/// ???????????? ???????????????????????? ?????????????????? ???????????? ???????????? ?????????????? ???????????? (????????????????????????)
		void CalcInfluenceFromVortexFarCells();


		void Gauss(std::vector<std::vector<double>>& A, const std::vector<PointsCopy>& pnt);

		/// ???????????? ?????????????? ???? ?????????? ???????????? ???????????? ???????????? ???? ?????????????? ???????? ?????????????? ?????????????? (?????? ?????????????? ????????)
		void CalcInfluenceFromPanels();
		void CalcInfluenceFromPanelsToPoints();
		/// ???????????????????? ?????????????? ???? ?????????? ???????????? ???????????? ???????????? ???? ?????????????? ???????? ?????????????? ?????????????? (?????? ?????????????? ????????)
		/// 
		/// ???????????????????? ???????????? ?????????????? IterativeInfluenceComputation ?????? ???????? ????????????????, ?????????? ???????????? 
		/// ???????????????????????? ???????? ???? ???????????????????????????? ????????????????, ?? ?????????????? ???? ???????????????????????? ?????????????? velSave
		void UpdateInfluence();

		void PushbackLowTrees();
	};

	//?????????????????? ?????????????????????? ??????????
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] });
	}

	// ?????????????????? a ???? ???????????????????? ????????????????????e ?? b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

}//nemaspace BH


#endif
