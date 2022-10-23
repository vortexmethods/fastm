/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.2    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/10/22     |
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
\version 1.2
\date 22 октября 2022 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_

#include <iostream>

namespace BH
{	
	/// Точность расчета скоростей:
	//  order = 1: монополь (0) + диполь (0)
	//  order = 2: монополь (0+1) + диполь (0) + квадруполь (0)
	//  order = 3: монополь (0+1+2) + диполь (0+1) + квадруполь (0) + октуполь (0)
	//  order = 4: монополь (0+1+2+3) + диполь (0+1+2) + квадруполь (0+1) + октуполь (0) + гексадекаполь (0)
	static const int order = 12;

	/// Максимальное количество уровней дерева
	static const int NumOfLevels = 10;

	/// Параметр точности 
	static const double theta = 1.25;

	/// Сравнение с прямым методом
	static const bool compare = true;

	/// Номер уровня дерева, до которого генерируются OMP нити
	static const int maxLevelOmp = -1;

	/// Число повторных запусков кода
	static const int runs = 3;

	//Признак подсчета числа операций
	//#define calcOp 


	//24
	const int codeLength = 14;
	const int twoPowCodeLength = (1 << codeLength);
	   
	static const double epsGMRES = 1e-15;

//#define CALCVORTEXVELO
#define CALCSHEET
//#define CALCVP


#ifdef CALCVORTEXVELO	
	#undef infToPanels
	#define needTreeVrt
	#undef needTreePan
	#undef needTreeVP

	// Имя файла со списком вихрей
	static const std::string vortexFile = "../../test/test10000.txt";

	// Радиус вихревого элемента
	static const double eps = 1e-3; //1e-4

	// Имя файла со списком вершин профиля
	static const std::string airfoilFile = "";	

	static const Point2D velInf = { 0.0, 0.0 };

	// Имя файла со списком точек вычисления скорости
	static const std::string vpFile = "";

	// Название задачи для файла с результатом
	static const std::string task = "10k";
#endif

#ifdef CALCSHEET	
	#define infToPanels
	
	#define needTreeVrt //для правой части
	#define needTreePan
	#undef needTreeVP

	#define linScheme // Признак линейной расчетной схемы 
	//#define asympScheme  // Признак асимптотической расчетной схемы
	
	// Имя файла со списком вихрей
	static const std::string vortexFile = "../../test/test1.txt";

	// Радиус вихревого элемента
	static const double eps = 0.032;
	//eps			N = 200;	N = 800;	N = 3200     N = 12800
	//Ellipse4x1:	0.0854014;	0.0214395;	0.00536141;  0.00134035;
	//Ellipse2x1:	0.0484234;	0.0121103;	0.00302764;
	//Circle:		2pi/n
	

	// Имя файла со списком вершин профиля
	static const std::string airfoilFile = "../../test/Wing/wing128.txt";
	//static const std::string airfoilFile = "../../test/Ellipse4x1/Ellipse12800.txt";
	
	// Скорость набегающего потока
	static const Point2D velInf = { 0.8660254037844386, -0.5 };

	// Имя файла со списком точек вычисления скорости
	static const std::string vpFile = "";

	// Название задачи для файла с результатом
	static const std::string task = "Wing128";
#endif

#ifdef CALCVP	
	#undef infToPanels

	#define needTreeVrt
	#define needTreePan
	#define needTreeVP

	#define linScheme // Признак линейной расчетной схемы 
	#define asympScheme  // Признак асимптотической расчетной схемы 

	// Имя файла со списком вихрей
	static const std::string vortexFile = "../../test/test0.txt";

	// Радиус вихревого элемента
	static const double eps = 0.032; //1e-4

	// Имя файла со списком вершин профиля
	static const std::string airfoilFile = "../../test/Wing/wing128.txt";
	//static const std::string airfoilFile = "../../test/Ellipse4x1/Ellipse12800.txt";

	// Скорость набегающего потока
	static const Point2D velInf = { 0.8660254037844386, -0.5 };

	// Имя файла со списком точек вычисления скорости
	static const std::string vpFile = "../../test/vp671.txt";

	// Название задачи для файла с результатом
	static const std::string task = "Wing128";
#endif

	
    

#ifdef asympScheme
	/*
	static const int nAngPoints = 2;
	static const int KK[] = { 0, 20 };
	static const int KKm[] = { 127, 19 };
	static const int p[] = { 1, 2 };
	static const int q[] = { 2, 6 };
	static const double mu[] = { 1.0 / 2.0, 1.0 / 3.0 };
	//*/

	//*
	static const int nAngPoints = 1;
	static const int KK[] = { 0 };
	static const int KKm[] = { 127 };
	static const int p[] = { 1 };
	static const int q[] = { 2 };
	static const double mu[] = { 1.0 / 2.0 };
	//*/

	/*
	static const int nAngPoints = 1;
	static const int KK[] = { 20 };
	static const int KKm[] = { 19 };
	static const int p[] = { 2 };
	static const int q[] = { 6 };
	static const double mu[] = { 1.0 / 3.0 };
	*/
#endif




	//КАК ПРАВИЛО, ПАРАМЕТРЫ НИЖЕ НЕ НУЖНО МЕНЯТЬ
	static const double PI = 3.1415926535897932384626;
	static const double IPI = 1.0 / PI;
	static const double IDPI = 0.5 / PI;

	// Радиус вихревого элемента (в квадрате)
	static const double eps2 = eps * eps;

	/// Минимальное количество вихрей в ячейке
	static const int minNumOfVort = 1;

	/// Минимальный размер ячейки
	static const double minDist = 0.0 * eps;

	// Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
	static const bool BSfromFile = true;

	// Сохранение скоростей по быстрому методу в файл
	static const bool save = true;




}//namespace BH

#endif