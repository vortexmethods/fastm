/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.10   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2021/05/17     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Point2D.h                                                        |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием классов Point2D, Point2Df
\author Марчевский Илья Константинович
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef POINT2D_H_
#define POINT2D_H_

//#ifndef __CUDACC__
//	#include "mpi.h"
//#endif

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<double, 2>, имеет дополнительные возможности:
	- поворота на заданный угол против часовой стрелки;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Point2D
		: public numvector<double, 2>
	{
	public:

//#ifndef __CUDACC__
//		/// MPI-описатель типа
//		static MPI_Datatype mpiPoint2D;
//#endif
		/// Пустой конструктор
		Point2D() { };

		/// \brief Конструктор и приведение типа из numvector<double, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 2>
		Point2D(const numvector<double, 2>& _r)
		{
			data()[0] = _r[0];
			data()[1] = _r[1];
		}//Point2D(...);

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		Point2D(const Point2D& _r)
		{
			data()[0] = _r[0];
			data()[1] = _r[1];
		}//Point2D(...)

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа double
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		Point2D(const std::initializer_list<double>& z)
		{
			for (size_t i = 0; i < 2; ++i)
				data()[i] = *(z.begin() + i);
		}//Point2D(...)
#endif

		/// Деструктор
		~Point2D() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
		/// \return новый вектор, полученный поворотом старого 
		Point2D rotated(const double angle = 1.5707963267948966192313216916398) const
		{
			Point2D res;
			double cosa = cos(angle);
			double sina = sin(angle);

			res[0] = data()[0] * cosa - data()[1] * sina;
			res[1] = data()[0] * sina + data()[1] * cosa;
			return res;
		}//rotated(...)

//#ifndef __CUDACC__
//		/// Cоздание MPI-описателя типа
//		static void CreateMpiType()
//		{
//			int          len[1] = { 2 };
//			MPI_Aint     pos[1] = { 0 };
//			MPI_Datatype typ[1] = { MPI_DOUBLE };
//
//			MPI_Type_create_struct(1, len, pos, typ, &mpiPoint2D);
//			MPI_Type_commit(&mpiPoint2D);
//
//
//			/*
//			int          len[3] = { 1, 2, 1 };
//			//MPI_Aint     pos[3] = { 0, offsetof(Point2D, data()), sizeof(Point2D) };
//			MPI_Aint     pos[3] = { 0, 0, sizeof(Point2D) };
//			MPI_Datatype typ[3] = { MPI_LB, MPI_DOUBLE, MPI_UB };
//
//			MPI_Type_create_struct(3, len, pos, typ, &mpiPoint2D);
//			MPI_Type_commit(&mpiPoint2D);
//			*/
//		}//CreateMpiType()
//#endif

		//operator numvector<double, 2>&() 
		//{
		//	return *this;
		//}
	};


	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<float, 2>, имеет дополнительные возможности:
	- поворота на заданный угол против часовой стрелки;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Point2Df
		: public numvector<float, 2>
	{
	public:

		//#ifndef __CUDACC__
		//		/// MPI-описатель типа
		//		static MPI_Datatype mpiPoint2Df;
		//#endif
				/// Пустой конструктор
		Point2Df() { };

		/// \brief Конструктор и приведение типа из numvector<float, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<float, 2>
		Point2Df(const numvector<float, 2>& _r)
		{
			data()[0] = _r[0];
			data()[1] = _r[1];
		}//Point2Df(...);

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		Point2Df(const Point2Df& _r)
		{
			data()[0] = _r[0];
			data()[1] = _r[1];
		}//Point2Df(...)

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа float
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		Point2Df(const std::initializer_list<float>& z)
		{
			for (size_t i = 0; i < 2; ++i)
				data()[i] = *(z.begin() + i);
		}//Point2Df(...)
#endif

		/// Деструктор
		~Point2Df() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
		/// \return новый вектор, полученный поворотом старого 
		Point2Df rotated(const float angle = 1.57079632f) const
		{
			Point2Df res;
			float cosa = cos(angle);
			float sina = sin(angle);

			res[0] = data()[0] * cosa - data()[1] * sina;
			res[1] = data()[0] * sina + data()[1] * cosa;
			return res;
		}//rotated(...)

//#ifndef __CUDACC__
//		/// Cоздание MPI-описателя типа
//		static void CreateMpiType()
//		{
//			int          len[1] = { 2 };
//			MPI_Aint     pos[1] = { 0 };
//			MPI_Datatype typ[1] = { MPI_FLOAT };
//
//			MPI_Type_create_struct(1, len, pos, typ, &mpiPoint2Df);
//			MPI_Type_commit(&mpiPoint2Df);
//
//
//			/*
//			int          len[3] = { 1, 2, 1 };
//			//MPI_Aint     pos[3] = { 0, offsetof(Point2Df, data()), sizeof(Point2Df) };
//			MPI_Aint     pos[3] = { 0, 0, sizeof(Point2Df) };
//			MPI_Datatype typ[3] = { MPI_LB, MPI_FLOAT, MPI_UB };
//
//			MPI_Type_create_struct(3, len, pos, typ, &mpiPoint2Df);
//			MPI_Type_commit(&mpiPoint2Df);
//			*/
//		}//CreateMpiType()
//#endif

		//operator numvector<float, 2>&() 
		//{
		//	return *this;
		//}
	};



}//namespace VMlib

using VMlib::Point2D;
using VMlib::Point2Df;

#endif
 