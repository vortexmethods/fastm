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
	/*!
	\brief Класс, определяющий основной алгоритм модификации метода Барнса - Хата
	\author Марчевский Илья Константинович
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\version 1.3
	\date 08 декабря 2022 г.
	*/
	class BarnesHut
	{
	public:
		///Ссылка на параметры, считываемые из файла
		const params& prm;
		
		///Список оберток положений вихрей
		std::vector<PointsCopy> pointsCopyVrt;
		
		///Список оберток положений центров панелей
		std::vector<PointsCopy> pointsCopyPan;
			
		///Список оберток положений точек вычисления скорости в области течения
		std::vector<PointsCopy> pointsCopyVP;
		
		///Умный yказатель на дерево вихрей
		mutable std::unique_ptr<MortonTree> treeVrt;
		
		///Умный yказатель на дерево центров панелей
		mutable std::unique_ptr<MortonTree> treePan;

		///Умный yказатель на дерево точек вычтсления скорости в области течения
		mutable std::unique_ptr<MortonTree> treeVP;

		
		/// \brief Конструктор для решения задачи NBODY о вычислении скоростей вихревых частиц
		///
		/// \param[in] prm константная ссылка на параметры, считываемые из файла
		/// \param[in] pointsVrt константная ссылка на список частиц
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt);
		
		/// \brief Конструктор для решения задачи BIE о решении граничного интегрального уравнения
		///
		/// \param[in] prm константная ссылка на параметры, считываемые из файла
		/// \param[in] pointsVrt константная ссылка на список частиц
		/// \param[in] pointsPan константная ссылка на список центров панелей
		/// \param[in] panPos константная ссылка на положения начал и концов панелей
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos);
		
		/// \brief Конструктор для решения задачи VP о вычислении скоростей в области течения
		///
		/// \param[in] prm константная ссылка на параметры, считываемые из файла
		/// \param[in] pointsVrt константная ссылка на список частиц
		/// \param[in] pointsPan константная ссылка на список центров панелей, в которых хранится постоянная составляющая решения на панелях
		/// \param[in] panPos константная ссылка на положения начал и концов панелей
		/// \param[in] sec константная ссылка на список, в котором хранятся линейные составляющие решения на панелях
		/// \param[in] pointsVP константная ссылка на список точек в области течения, где рассчитываются скорости
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt, const std::vector<Vortex2D>& pointsPan, const std::vector<Point2D>& panPos, const std::vector<double>& sec, const std::vector<Vortex2D>& pointsVP);
		
		/// \brief Формирование элемента списка pointsCopyPan, соответствующего панели
		///
		/// \param[in] panCenter константная ссылка на центр панели, в котором хранится постоянная составляющая решения на панели
		/// \param[in] panBegin константная ссылка на положение начала панели
		/// \param[in] panEnd константная ссылка на положение конца панели
		/// \param[in] gamLin линейная составляющая решения на панели
		void CreatePan(const Vortex2D& panCenter, const Point2D& panBegin, const Point2D& panEnd, double gamLin);
		
		/// \brief Преобразование панели к панели, на которой задана асимптотика решения (в окрестности угловой точки)
		///
		/// \param[in,out] panel ссылка на преобразовываемую панель
		/// \param[in] mu показатель степени в особенности решения в угловой точке
		/// \param[in] infinityAtBegin значение true, если особенность в начале панели, и false, если особенность в конце панели
		void ConvertToAsypmPanel(PointsCopy& panel, double mu, bool infinityAtBegin);

		/// Деструктор
		~BarnesHut() {};

		/// \brief Построение одного дерева tree на основе заданных точек pointsCopy  
		///
		/// \param[in,out] tree умный указатель на дерево
		/// \param[in] maxTreeLevel максимальная глубина при обходе дерева
		/// \param[in] pointsCopy неконстантная ссылка на список данных (оберток) точек, по которым строится дерево
		/// \param[in] ifpan признак того, что дерево состоит из центров панелей, а не из отдельных частиц
		/// \param[in,out] time время, затрачиваемое на построение дерева (накопительный итог)
		void BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, bool ifpan, double& time);

		/// \brief Построение всех нужных деревьев на основе заданных точек pointsCopy  
		///
		/// \param[in,out] time время, затрачиваемое на построение дерева (накопительный итог)
		void BuildNecessaryTrees(double& time);

		/// \brief Расчет влияния в точках дерева, характерных для решаемой задачи (определяется внутри функции) 
		/// 
		/// \param[out] result ссылка на вектор, в который сохраняются вычисленные скорости
		/// \param[in,out] timeParams время расчета параметров деревьев (накопительный итог)
		///	\param[in,out] timeInfl время расчета влияния (накопительный итог)
		void InfluenceComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl);
		
		/// \brief Расчет правой части при решении ГИУ
		/// 
		/// \param[out] result ссылка на вектор, в который сохраняются вычисленные скорости
		/// \param[in,out] timeParams время расчета параметров дерева (накопительный итог)
		///	\param[in,out] timeInfl время расчета влияния (накопительный итог)
		void RhsComputation(std::vector<Point2D>& result, double& timeParams, double& timeInfl);
	
		/// \brief Расчет влияния на панели от них же внутри итерационного алгоритма
		/// 
		/// \param[out] result ссылка на вектор, в который сохраняются вычисленные скорости
		/// \param[in] newGam константная ссылка на вектор текущих "циркуляций" на панелях
		/// \param[in,out] timeParams время расчета параметров дерева (накопительный итог)
		///	\param[in,out] timeInfl время расчета влияния (накопительный итог)
		void IterativeInfluenceComputation(std::vector<Point2D>& result, const std::vector<double>& newGam, double& timeParams, double& timeInfl);

		/// \brief Расчет влияния вихревого следа на панели при решении ГИУ 
		/// 
		/// \param[out] rhs ссылка на заполняемый вектор правой части
		/// \param[in,out] timeParams время расчета параметров дерева (накопительный итог)
		///	\param[in,out] timeInfl время расчета влияния (накопительный итог)
		void FillRhs(std::vector<double>& rhs, double& timeParams, double& timeInfl);

		/// \brief Обновление циркуляций вихревых элементов (для итерационного решения СЛАУ/ГИУ)
		/// 
		/// \param[in] newGam константная ссылка на вектор нового приближения (поправки, невязки, ...) 
		void UpdateGams(const std::vector<double>& newGam);

	};

}//namespace BH

#endif
