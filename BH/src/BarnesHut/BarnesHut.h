/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
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
\version 1.3
\date 08 декабря 2022 г.
*/


#ifndef BARNESHUT_H_
#define BARNESHUT_H_

#include "Tree.h"

namespace BH
{

	class BarnesHut
	{
	public:
		///Ссылка на параметры
		const params& prm;
		
		///Список точек: для вихрей и для панелей
		std::vector<PointsCopy> pointsCopyPan, pointsCopyVrt, pointsCopyVP;
		
		/// умный yказатель на дерево 
		mutable std::unique_ptr<MortonTree> treePan, treeVrt, treeVP;

		/// Конструктор	
		///
		/// \param[in] pointsPan константная ссылка на список панелей
		/// \param[in] pointsVrt константная ссылка на список частиц
		/// \param[in] panPos константная ссылка на начала и концы панелей (для решения СЛАУ)
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos);
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos, std::vector<double> sec, const std::vector<Vortex2D>& pointsVP);
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt);

		void CreatePan(const Vortex2D& panCenter, const Point2D& panBegin, const Point2D& panEnd, double gamLin);
		void ConvertToAsypmPanel(PointsCopy& panel, double mu, bool infinityAtBegin);

		/// Деструктор
		~BarnesHut() {};

		///\brief Построение одного дерева tree на основе заданных точек pointsCopy  
		void BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, bool ifpan, double& time);

		///\brief Построение всех нужных деревьев на основе заданных точек pointsCopy  
		void BuildNecessaryTrees(double& time);

		///\brief Расчет влияния
		/// \param[out] result --- вектор из скоростей точек для сохранения результата
		/// \param[out] timeParams --- время расчета параметров дерева
		///	\param[out] timeInfl --- время расчета влияния
		void InfluenceComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl);
		void RhsComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl);
	
		/// \brief Расчет влияния точек result самих на себя внутри итерационного алгоритма
		/// \param[out] result --- вектор из скоростей точек
		/// \param[in] newGam --- вектор нового приближения/поправки/невязки
		void IterativeInfluenceComputation(std::vector<Point2D>& result, std::vector<double>& newGam, double& timeParams, double& timeInfl);

		/// \brief Расчет влияния вихревого следа на правую часть СЛАУ 
		/// \param[out] rhs --- вектор правой части
		void FillRhs(std::vector<double>& rhs, double& timeParams, double& timeInfl);

		/// Обновление циркуляций вихревых элементов (для решения СЛАУ)
		/// \param[in] newGam --- вектор нового приближения/поправки/невязки 
		void UpdateGams(std::vector<double>& newGam);

	};

}//namespace BH

#endif
