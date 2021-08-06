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
\brief Заголовок класса Tree
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
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
		/// Габариты прямоугольника
		double iLeft, iRight, iTop, iBottom;

		/// Центр прямоугольника
		Point2D posCentre;

		/// Номер моего уровня (глубина дерева)
		int level;
		/// Номер потомка внутри уровня
		int p;

		/// Являюсь ли я нижним уровнем?
		bool lowLevel = false;

		/// Потомки данного уровня дерева
		std::unique_ptr<Tree> mChildren[2];

		/// Монопольный момент ячейки 
		double mon;
		/// Дипольный, квадрупольный, октупольный и гексадекапольный моменты ячейки 
		Point2D dip, qua, oct, hex;

		/// \brief Коэффициенты для вычисления скоростей
		numvector<Point2D, 4> E;

	public:

		std::vector<PointsCopy>& points;

		/// Вектор указателей на деревья в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для нижних уровней
		std::vector<Tree*> closeTrees;

		/// Вектор указателей на ячейки нижнего уровня 
		//имеет смысл только для корня
		std::vector<Tree*> lowTrees;

		/// Индексы элементов в исходном массиве точек
		std::vector<int> index;

		/// Указатель на корень (для заполнения lowTrees)
		Tree* BigTree;

		std::vector<std::unique_ptr<Tree>> treePtr;

		/// Конструктор
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

		/// Деструктор
		~Tree() {};


		/// Построение корня дерева на основе заданных вихрей
		///
		/// \param[in] points ссылка на список вихрей, на основе которых строится прямоугольник
		void makeTree(std::vector<PointsCopy>& points);

		/// Построение структуры дерева
		void CreateTree();

		/// Вычисление параметров дерева (циркуляций и центров завихренности)	
		// для нижнего уровня считаем по набору вихрей, для остальных --- по потомкам
		void CalculateTreeParams();
		/// Вычисление параметров нижней ячейки по набору вихрей (для расчета скоростей)
		void CalcPointsParams();
		/// Вычисление параметров нижней ячейки по набору панелей (для решения системы)
		void CalcPanelsParams();

		/// Расчет конвективных скоростей вихрей внутри данной ячейки нижнего уровня от вихрей всех ближних уровней
		void CalcConvVeloByBiotSavart();

		/// Обнуление коэффициентов локальных разложений
		void ClearCoeff();

		/// Расчет коэффициентов разложения в ряд Тейлора внутри ячейки нижнего уровня
		///
		/// \param[in] lowTree указатель на дерево нижнего уровеня
		/// \param[in] calcCloseTrees = true (по умолчанию), если нужно заполнить вектор ближних ячеек
		void CalcLocalCoeffToLowLevel(Tree* lowTree, bool calcCloseTrees = true);

		/// Расчет конвективных скоростей внутри дерева нижнего уровня (приближенный)
		void CalcInfluenceFromVortexFarCells();


		void Gauss(std::vector<std::vector<double>>& A, const std::vector<PointsCopy>& pnt);

		/// Расчет влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
		void CalcInfluenceFromPanels();
		void CalcInfluenceFromPanelsToPoints();
		/// Обновление влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
		/// 
		/// Вызывается внутри функции IterativeInfluenceComputation для всех итераций, кроме первой 
		/// Коэффициенты СЛАУ не рассчитываются напрямую, а берутся из сохраненного массива velSave
		void UpdateInfluence();

		void PushbackLowTrees();
	};

	//умножение комплексных чисел
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] });
	}

	// умножение a на комплексно сопряженноe к b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

}//nemaspace BH


#endif
