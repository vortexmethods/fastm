/*--------------------------------*- BHcu -*-----------------*---------------*\
| #####   ##  ##                |                            | Version 1.1    |
| ##  ##  ##  ##   ####  ##  ## |  BHcu: Barnes-Hut method   | 2022/08/28     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
| Info: Source code of BHcu                                                   |
|                                                                             |
| This file is part of BHcu.                                                  |
| BHcu is free software: you can redistribute it and/or modify it             |
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
| along with BHcu.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Параметры решаемой задачи
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 28 августа 2022 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_

#include <string>

namespace BHcu
{
	/// Радиус вихревого элемента
#define EPS 1e-4
/// Параметр точности 
#define THETA 1.00

	static const int order = 11;   //1-MONOPOLE
	
	// Имя файла с задачей
	static const std::string nameFile = "../../test/Wake/wake2m.txt";
	// Название задачи для файла с результатом
	static const std::string task = "2m";

	//Сравнение с прямым методом
	static const bool compare = true;


	//Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
	static const bool BSfromFile = true;
	//Сохранение скоростей по быстрому методу в файл
	static const bool save = false;

	//Число повторных запусков кода
	static const int runs = 10;

	//Тип данных 
#define CALCinDOUBLE

//Печатать ли полную информацию о device
	static const bool printFullCUDAinfo = false;



	//КАК ПРАВИЛО, ПАРАМЕТРЫ НИЖЕ НЕ НУЖНО МЕНЯТЬ

	static const int dev = 0;

#ifdef CALCinFLOAT
#define real float
#define real2 float2
#define realmax fmaxf
#define realmin fminf
#define realPoint Point2Df
#define realVortex Vortex2Df
#endif

#ifdef CALCinDOUBLE
#define real double
#define real2 double2
#define realmax fmax
#define realmin fmin
#define realPoint Point2D
#define realVortex Vortex2D
#endif

}//namespace BHcu

#endif