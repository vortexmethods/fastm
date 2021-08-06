/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
| Info: Source code of FFT                                                    |
|                                                                             |
| This file is part of FFT.                                                   |
| FFT is free software: you can redistribute it and/or modify it              |
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
| along with FFT.  If not, see <http://www.gnu.org/licenses/>.                |
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

namespace FFT
{

#include <iostream>
	// радиус вихревого элемента
	static const double eps = 1e-3;

	// Имя файла с задачей
	static const std::string nameFile = "../../test/test100000.txt";
	// Название задачи для файла с результатом
	static const std::string task = "100k";

	// Число узлов сетки
	static const int NP = 128;
	// Количество окаймляющих слоев в ближней зоне
	static const int sizeOfCloseZone = 3;
	// Тип интерполяции: 0 - билинейная, 1 - эрмитова
	static const int interpType = 1;

	// Количество OMP нитей для 1 и 3 операций соответственно (0 - максимально возможное)
	static const int numTh1 = 4;
	static const int numTh3 = 4;

	//Сравнение с прямым методом
	static const bool compare = true;

	//Число повторных запусков кода
	static const int runs = 10;


	//КАК ПРАВИЛО, ПАРАМЕТРЫ НИЖЕ НЕ НУЖНО МЕНЯТЬ

	// радиус вихревого элемента (в квадрате)
	static const double eps2 = eps * eps;

	//Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
	static const bool BSfromFile = true;
	//Сохранение скоростей по быстрому методу в файл
	static const bool save = false;

}//namespace FFT

#endif