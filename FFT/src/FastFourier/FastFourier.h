/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: FastFourier.h                                                    |
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
\brief Описание основного класса
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
*/

#ifndef FASTFOURIER_H_
#define FASTFOURIER_H_

#include <vector>

//#define EIGEN_FFTW_DEFAULT
#include "Eigen/Eigen"
#include "unsupported/Eigen/FFT"

#include "Params.h"
#include "Vortex2D.h"


namespace FFT
{

	const double PI = 3.1415926535897932;
	const double IDPI = 0.15915494309189534;

	/// Описание структуры одной ячейки сетки
	struct cell
	{
		/// Координаты левого нижнего узла ячейки
		Point2D coord;

		/// Вектор из указателей на ячейки ближней зоны
		std::vector<cell*> closeZone;

		/// Вектор номеров реальных ячеек ближней зоны в общем массиве
		std::vector<int> numOfCloseZone;

		/// Вектор коррекционных скоростей в узлах ячейки 
		numvector<Point2D, 4> velCorr;

		/// Вектор производных скорости в узлах ячейки 
		numvector<Point2D, 4> velDerX, velDerY;

		/// Индексы вихрей данной ячейки в общем массиве
		std::vector<int> index;

		/// Индексы вихрей, попадающих в "теневую" ближнюю зону (справа и слева)
		std::vector<int> indexRight, indexLeft;

		/// Вектор циркуляций в 9-ти соседних ячейках, на которые "размазывается" влияние данной
		numvector<double, 16> gam16;

	};



	class FastFourier
	{

	private:
		/// Длина "ближней зоны" (в ячейках)
		int lengthOfCloseZone;

		/// Вектор координат едининых вихрей, используемых для построения матрицы коррекции
		std::vector<Point2D> posForC;

		/// Матрица коррекции С
		//std::vector<numvector<double, 16>> Cmatr; 
		std::vector<double> Cmatr;

		/// Число узлов сетки (Nx - по горизонтали, Ny - по вертикали )
		const numvector<int, 2> nNodes;

		/// Шаги сетки
		const Point2D h;

		/// Координаты левого нижнего угла области (const)
		const Point2D lD;

		/// Координаты правого верхнего угла области (const)
		Point2D rU;

		/// Размеры области
		// dim[0] --- по горизонтали 
		// dim[1] --- по вертикали
		const Point2D dim;

		/// Тип интерполяции с узлов сетки на ВЭ (0 --- билинейная, 1 --- интерполяция Эрмита)
		const int interpType;

	public:
		/// Вектор из вихрей
		std::vector<Vortex2D>& points;

		/// Вектор из скоростей
		std::vector<Point2D>& velo;

		/// Вектор из ячеек сетки
		std::vector<cell> cells;

		/// Вектор циркуляций в узлах сетки (свой для каждого процесса MPI)
		std::vector<double> gam;

		/// Вектор скоростей, полученных решением уравнения Пуассона на всей сетке 
		std::vector<Point2D> velPuasson;

		std::vector<double> bkx;
		std::vector<double> bky;
		std::vector<double> omega;
		std::vector<double> hUx;
		std::vector<double> hUy;


		/// Конструктор инициализации
		FastFourier(std::vector<Vortex2D>& wake, std::vector<Point2D>& _velo, Point2D _lD, Point2D _dim, numvector<int, 2> _nNodes, int _interpType);

		/// Одномерный оператор сглаживания 
		double M4(double x);

		/// Двумерный оператор сглаживания
		double W(const Point2D& r);

		/// Вычисление вектора циркуляций в узлах сетки с помощью оператора сглаживания
		void FillGam();

		/// Инициализация ячеек сетки
		void InitializationCells();

		/// Вычисление матрицы коррекции С 
		// Считывание из файла
		void ReadMatrC();

		/// Вычисление коррекционных скоростей и производных, заполнение Velcorr, VelDerX и VelDerY для каждой ячейки
		void CalcCorrectionVel();

		/// Расчет скоростей ВЭ внутри q ячейки по закону Био --- Савара от ближней зоны
		inline void CalcVeloByBiotSavartForCell(const cell& q);

		/// Функция ядра Био --- Савара
		Point2D Skos(const Point2D& R, const Point2D& X) const;

		/// Решение уравнение Пуассона, заполнение вектора VelPuasson
		void PuassonSolver(double& tExclude);

		/// Линейная интерполяция с узлов на ВЭ в q-й ячейке
		void InterpolationFromCell(cell& q);

		/// Интерполяция Эрмита с узлов на ВЭ в q-й ячейке
		void HermiteInterpolationFromCell(cell& q);

		/// Расчет итоговых скоростей
		void CalcTotalVelo();

		/// Расчет циркуляций в 16 ближайших узлах с помощью оператора Монагана
		numvector<double, 16> CalcGam16(const cell& q);

		/// Деструктор
		~FastFourier() {};
	};

}//namespace FFT

#endif