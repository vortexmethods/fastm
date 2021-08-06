/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.10   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2021/05/17     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Vortex2D.h                                                       |
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
\brief Заголовочный файл с описанием классов Vortex2D, Vortex2Df
\author Марчевский Илья Константинович
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef VORTEX2D_H_
#define VORTEX2D_H_

#include "Point2D.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий двумерный вихревой элемент	
	\author Марчевский Илья Константинович
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Vortex2D
	{
	private:
		/// Радиус-вектор вихря
		Point2D pos;

		/// Циркуляция вихря
		double gam;

	public:

//#ifndef __CUDACC__
//		/// MPI-описатель типа
//		static MPI_Datatype mpiVortex2D;
//#endif
//
//		static size_t offsPos;
//		static size_t offsGam;

		/// Пустой конструктор
		Vortex2D() {};

		/// \brief Конструктор инициализации
		///
		/// \param[in] _r константная ссылка на радиус-вектор положения вихря
		/// \param[in] _g циркуляция (интенсивность) вихря
		Vortex2D(const Point2D& _r, const double _g)
			: pos(_r), gam(_g) { }

		/// Деструктор
		~Vortex2D() {};

		/// \brief Функция для доступа к радиус-вектору вихря
		/// \return ссылка на радиус-вектор вихря
		Point2D& r() { return pos; }

		/// \brief Функция для доступа для чтения к радиус-вектору вихря
		/// \return константная ссылка на радиус-вектор вихря
		const Point2D& r() const { return pos; }

		/// \brief Функция для доступа к циркуляции вихря
		/// \return ссылка на циркуляцию вихря
		double& g() { return gam; }

		/// \brief Функция для доступа для чтения к циркуляции вихря
		/// \return константная ссылка на циркуляцию вихря
		const double& g() const { return gam; }

//#ifndef __CUDACC__
//		/// Cоздание MPI-описателя типа
//		static void CreateMpiType()
//		{
//			int          len[2] = { 2, 1 };
//			MPI_Aint     pos[2] = { offsetof(Vortex2D, pos), offsetof(Vortex2D, gam) };
//			MPI_Datatype typ[2] = { MPI_DOUBLE, MPI_DOUBLE };
//
//			MPI_Type_create_struct(2, len, pos, typ, &mpiVortex2D);
//			MPI_Type_commit(&mpiVortex2D);
//
//			/*
//			int          len[4] = { 1, 2, 1, 1 };
//			MPI_Aint     pos[4] = { 0, offsetof(Vortex2D, pos), offsetof(Vortex2D, gam), sizeof(Vortex2D) };
//			MPI_Datatype typ[4] = { MPI_LB, MPI_DOUBLE, MPI_DOUBLE, MPI_UB };
//
//			MPI_Type_create_struct(4, len, pos, typ, &mpiVortex2D);
//			MPI_Type_commit(&mpiVortex2D);
//			*/
//
//			offsPos = offsetof(Vortex2D, pos);
//			offsGam = offsetof(Vortex2D, gam);
//		}//CreateMpiType()
//#endif
	};

	/// Определение типа данных - источника, имеющего ту же структуру, что и вихрь
	typedef Vortex2D Source2D;





	/*!
	\brief Класс, опеделяющий двумерный вихревой элемент (float)
	\author Марчевский Илья Константинович
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Vortex2Df
	{
	private:
		/// Радиус-вектор вихря
		Point2Df pos;

		/// Циркуляция вихря
		float gam;

	public:

		//#ifndef __CUDACC__
		//		/// MPI-описатель типа
		//		static MPI_Datatype mpiVortex2Df;
		//#endif
		//
		//		static size_t offsPos;
		//		static size_t offsGam;

				/// Пустой конструктор
		Vortex2Df() {};

		/// \brief Конструктор инициализации
		///
		/// \param[in] _r константная ссылка на радиус-вектор положения вихря
		/// \param[in] _g циркуляция (интенсивность) вихря
		Vortex2Df(const Point2Df& _r, const float _g)
			: pos(_r), gam(_g) { }

		/// Деструктор
		~Vortex2Df() {};

		/// \brief Функция для доступа к радиус-вектору вихря
		/// \return ссылка на радиус-вектор вихря
		Point2Df& r() { return pos; }

		/// \brief Функция для доступа для чтения к радиус-вектору вихря
		/// \return константная ссылка на радиус-вектор вихря
		const Point2Df& r() const { return pos; }

		/// \brief Функция для доступа к циркуляции вихря
		/// \return ссылка на циркуляцию вихря
		float& g() { return gam; }

		/// \brief Функция для доступа для чтения к циркуляции вихря
		/// \return константная ссылка на циркуляцию вихря
		const float& g() const { return gam; }

		//#ifndef __CUDACC__
		//		/// Cоздание MPI-описателя типа
		//		static void CreateMpiType()
		//		{
		//			int          len[2] = { 2, 1 };
		//			MPI_Aint     pos[2] = { offsetof(Vortex2Df, pos), offsetof(Vortex2Df, gam) };
		//			MPI_Datatype typ[2] = { MPI_FLOAT, MPI_FLOAT };
		//
		//			MPI_Type_create_struct(2, len, pos, typ, &mpiVortex2Df);
		//			MPI_Type_commit(&mpiVortex2Df);
		//
		//			/*
		//			int          len[4] = { 1, 2, 1, 1 };
		//			MPI_Aint     pos[4] = { 0, offsetof(Vortex2Df, pos), offsetof(Vortex2Df, gam), sizeof(Vortex2Df) };
		//			MPI_Datatype typ[4] = { MPI_LB, MPI_FLOAT, MPI_FLOAT, MPI_UB };
		//
		//			MPI_Type_create_struct(4, len, pos, typ, &mpiVortex2Df);
		//			MPI_Type_commit(&mpiVortex2Df);
		//			*/
		//
		//			offsPos = offsetof(Vortex2Df, pos);
		//			offsGam = offsetof(Vortex2Df, gam);
		//		}//CreateMpiType()
		//#endif
	};

	/// Определение типа данных - источника, имеющего ту же структуру, что и вихрь
	typedef Vortex2D Source2D;
	typedef Vortex2Df Source2Df;


}//namespace VMlib

using VMlib::Vortex2D;
using VMlib::Vortex2Df;
#endif
 
