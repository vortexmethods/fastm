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
\brief Заголовок класса-обертки для точек и панелей
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
/*!
\brief Класс-обертка для точек и панелей
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

	class PointsCopy : public Vortex2D
	{
	public:

		/// Вычисленная скорость
		Point2D veloCopy; 
		
		/// Скорость, спроецированная на линейную базисную функцию (для кусочно-линейной схемы)
		Point2D veloCopyLin;

		/// Проекция вычисленной скорости на касательную к профилю
		double velTau;
		
		/// Проекция скорости, спроецированная на линейную базисную функцию (для кусочно-линейной схемы), на касательную к профилю
		double velTauLin;

		/// Компоненты матрицы скосов
		std::vector<Point2D> i00;

		/// Компоненты предобуславливателя в кусочно-постоянной схеме T0
		Point2D a, c;

#ifdef linScheme
		/// Компоненты матрицы скосов
		std::vector<Point2D> i01, i10, i11;
		
		/// Компоненты предобуславливателя в кусочно - постоянной схеме T0
		Point2D a1, c1; //для предобуславливателя в Т1
		
		/// Линейная составляющая решения на панели (к-т при линейной базисной функции)
		double gamLin;
#endif

		/// Начало панели
		Point2D panBegin;
		
		/// Конец панели
		Point2D panEnd;

		/// Мультипольные моменты 
		std::vector<Point2D> Mpan;

#ifdef asympScheme
		/// Мультипольные моменты в схеме с выделением асимптотики решения
		std::vector<Point2D> MpanAs;
#endif

		/// Орт касательной к панели
		Point2D tau;

		/// Длина панели
		double len;

		/// Инициализирующий конструктор
		PointsCopy() :
			Vortex2D({ 0.0, 0.0 }, 0.0),
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 })
		{};

		/// Конструктор копирования из точки
		PointsCopy(const Point2D& r) :
			Vortex2D(r, 0.0),
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 })
		{};

		/// Конструктор копирования из вихря
		PointsCopy(const Vortex2D& vtx_) :
			Vortex2D(vtx_), veloCopy({ 0.0, 0.0 }), veloCopyLin({ 0.0, 0.0 }), 
			a({ 0.0, 0.0 }), c({ 0.0, 0.0 })
		{};

		/// Деструктор
		~PointsCopy() {};
	};

}//namespace BH

#endif
