/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
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
\brief Параметры решаемой задачи
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

#ifndef DEF_H_
#define DEF_H_

#include <string>

namespace FMM
{

	//Число членов разложения
	const int nt = 10;

	// радиус вихревого элемента
	static const double eps = 1e-3;

	// Имя файла с задачей
	static const std::string nameFile = "../../test/test2000000.txt";
	// Название задачи для файла с результатом
	static const std::string task = "2m";

	// Число уровней дерева 
	const int maxLevel = 8;

	//Сравнение с прямым методом
	static const bool compare = true;

	//Число повторных запусков кода
	static const int runs = 5;

	//Признак подсчета числа операций
	//#define calcOp


	//КАК ПРАВИЛО, ПАРАМЕТРЫ НИЖЕ НЕ НУЖНО МЕНЯТЬ

	// радиус вихревого элемента (в квадрате)
	static const double eps2 = eps * eps;

	// вычисление потенциала
	const bool calcPot = false;

	//Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
	static const bool BSfromFile = true;

	//Сохранение скоростей по быстрому методу в файл
	static const bool save = false;

	static int n; 

}//namespace FMM

#endif
