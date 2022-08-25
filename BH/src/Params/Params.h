/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
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
\brief Параметры решаемой задачи
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_

#include <iostream>

namespace BH
{

	// радиус вихревого элемента
	static const double eps = 1e-3;

	// Имя файла с задачей
	static const std::string nameFile = "../../test/test2000000.txt";
	// Название задачи для файла с результатом
	static const std::string task = "2m";

	/// Точность расчета скоростей:
	//  для решения ГИУ доступны order = 1-4, для скоростей order = 1-6
	//  order = 1: монополь (0) + диполь (0)
	//  order = 2: монополь (0+1) + диполь (0) + квадруполь (0)
	//  order = 3: монополь (0+1+2) + диполь (0+1) + квадруполь (0) + октуполь (0)
	//  order = 4: монополь (0+1+2+3) + диполь (0+1+2) + квадруполь (0+1) + октуполь (0) + гексадекаполь (0)
	static const int order = 6;
 // 1 - dipole

	/// Максимальное количество уровней дерева
	static const int NumOfLevels = 17;

	/// Параметр точности 
	static const double theta = 1.0;

	/// Сравнение с прямым методом
	static const bool compare = true;

	/// Номер уровня дерева, до которого генерируются OMP нити
	static const int maxLevelOmp = 4;

	/// Число повторных запусков кода
	static const int runs = 10;

	//Признак подсчета числа операций
	//#define calcOp 


	//24
	const int codeLength = 10;
	const int twoPowCodeLength = (1 << codeLength);


	//КАК ПРАВИЛО, ПАРАМЕТРЫ НИЖЕ НЕ НУЖНО МЕНЯТЬ

	// радиус вихревого элемента (в квадрате)
	static const double eps2 = eps * eps;

	/// Минимальное количество вихрей в ячейке
	static const int minNumOfVort = 1;

	/// Минимальный размер ячейки
	static const double minDist = 0.0 * eps;

	//Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
	static const bool BSfromFile = true;

	//Сохранение скоростей по быстрому методу в файл
	static const bool save = false;

	//#define pan 	
	//#define linScheme

	static const double PI = 3.1415926535897932384626;
	static const double IPI = 1.0 / PI;
	static const double IDPI = 0.5 / PI;

}//namespace BH

#endif