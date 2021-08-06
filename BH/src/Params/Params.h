/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.0    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2021/08/05     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
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
\version 1.0
\date 05 августа 2021 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_

#include <iostream>

namespace BH
{

	// радиус вихревого элемента
	static const double eps = 1e-3;

	// Имя файла с задачей
	static const std::string nameFile = "../../test/test100000.txt";
	// Название задачи для файла с результатом
	static const std::string task = "100k";

	///  Точность расчета скоростей:
	// для решения ГИУ доступны order = 1-4, для скоростей order = 1-6
	//  order = 1: монополь (0) + диполь (0)
	//  order = 2: монополь (0+1) + диполь (0) + квадруполь (0)
	//  order = 3: монополь (0+1+2) + диполь (0+1) + квадруполь (0) + октуполь (0)
	//  order = 4: монополь (0+1+2+3) + диполь (0+1+2) + квадруполь (0+1) + октуполь (0) + гексадекаполь (0)
	static const int order = 4;

	/// Максимальное количество уровней дерева
	static const int NumOfLevels = 14;

	/// Параметр точности 
	static const double theta = 0.5;

	//Сравнение с прямым методом
	static const bool compare = true;

	/// Номер уровня дерева, до которого генерируются OMP нити
	static const int maxLevelOmp = 1;

	//Число повторных запусков кода
	static const int runs = 10;

	//Признак подсчета числа операций
	//#define calcOp 


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

}//namespace BH

#endif