/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: BarnesHut.h                                                      |
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
\brief Заголовок основного класса BarnesHut
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/


#ifndef BARNESHUT_H_
#define BARNESHUT_H_

#include <iostream>
#include "Tree.h"

namespace BH
{

	class BarnesHut
	{
	public:

		double tStart, tFin;
		std::vector<PointsCopy> pointsCopy;
		
		/// умный yказатель на дерево 
		mutable std::unique_ptr<MortonTree> tree;

		/// Конструктор	
		///
		/// \param[in] points константная ссылка на список элементов
		/// \param[in] panPos константная ссылка на начала и концы панелей (для решения СЛАУ)
		BarnesHut(const std::vector<Vortex2D>& points, const std::vector<Point2D>& panPos);
		BarnesHut(const std::vector<Vortex2D>& points);

		/// Деструктор
		~BarnesHut() {};

		void BuildTree(double& time);

		///\brief Расчет влияния точек result самих на себя
		/// \param[out] result --- вектор из скоростей точек 
		/// \param[in] type --- тип расчета влияния "напрямую"
		// 0 - Био --- Савар (по умолчанию), схема точка-точка БЕЗ РЕЗЕРВА
		// 1 - Схема панель-панель
		// 2 - Схема панель-точка
		void InfluenceComputation(std::vector<Point2D>& result, int type, double& timeParams, double& timeInfl);

		/// \brief Расчет влияния точек result самих на себя внутри итерационного алгоритма
		/// \param[out] result --- вектор из скоростей точек
		/// \param[in] newGam --- вектор нового приближения/поправки/невязки
		/// \param[in] type --- тип расчета влияния "напрямую"
		// 0 - Био --- Савар (по умолчанию), схема точка-точка
		// 1 - Схема панель-панель
		// 2 - Схема панель-точка
		
		//void IterativeInfluenceComputation(std::vector<Point2D>& result, std::vector<double>& newGam, int type = 0);

		/// Обновление циркуляций вихревых элементов (для решения СЛАУ)
		void UpdateGams(std::vector<double>& newGam);

		/// Обнуление временной статистики
		void ClearTimestat();

	};

}//namespace BH

#endif
