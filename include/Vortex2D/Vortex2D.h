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
\brief ???????????????????????? ???????? ?? ?????????????????? ?????????????? Vortex2D, Vortex2Df
\author ???????????????????? ???????? ????????????????????????????
\version 1.10
\date 17 ?????? 2021 ??.
*/

#ifndef VORTEX2D_H_
#define VORTEX2D_H_

#include "Point2D.h"

namespace VMlib
{

	/*!
	\brief ??????????, ?????????????????????? ?????????????????? ???????????????? ??????????????	
	\author ???????????????????? ???????? ????????????????????????????
	\version 1.10
	\date 17 ?????? 2021 ??.
	*/
	class Vortex2D
	{
	private:
		/// ????????????-???????????? ??????????
		Point2D pos;

		/// ???????????????????? ??????????
		double gam;

	public:

//#ifndef __CUDACC__
//		/// MPI-?????????????????? ????????
//		static MPI_Datatype mpiVortex2D;
//#endif
//
//		static size_t offsPos;
//		static size_t offsGam;

		/// ???????????? ??????????????????????
		Vortex2D() {};

		/// \brief ?????????????????????? ??????????????????????????
		///
		/// \param[in] _r ?????????????????????? ???????????? ???? ????????????-???????????? ?????????????????? ??????????
		/// \param[in] _g ???????????????????? (??????????????????????????) ??????????
		Vortex2D(const Point2D& _r, const double _g)
			: pos(_r), gam(_g) { }

		/// ????????????????????
		~Vortex2D() {};

		/// \brief ?????????????? ?????? ?????????????? ?? ????????????-?????????????? ??????????
		/// \return ???????????? ???? ????????????-???????????? ??????????
		Point2D& r() { return pos; }

		/// \brief ?????????????? ?????? ?????????????? ?????? ???????????? ?? ????????????-?????????????? ??????????
		/// \return ?????????????????????? ???????????? ???? ????????????-???????????? ??????????
		const Point2D& r() const { return pos; }

		/// \brief ?????????????? ?????? ?????????????? ?? ???????????????????? ??????????
		/// \return ???????????? ???? ???????????????????? ??????????
		double& g() { return gam; }

		/// \brief ?????????????? ?????? ?????????????? ?????? ???????????? ?? ???????????????????? ??????????
		/// \return ?????????????????????? ???????????? ???? ???????????????????? ??????????
		const double& g() const { return gam; }

//#ifndef __CUDACC__
//		/// C?????????????? MPI-?????????????????? ????????
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

	/// ?????????????????????? ???????? ???????????? - ??????????????????, ???????????????? ???? ???? ??????????????????, ?????? ?? ??????????
	typedef Vortex2D Source2D;





	/*!
	\brief ??????????, ?????????????????????? ?????????????????? ???????????????? ?????????????? (float)
	\author ???????????????????? ???????? ????????????????????????????
	\version 1.10
	\date 17 ?????? 2021 ??.
	*/
	class Vortex2Df
	{
	private:
		/// ????????????-???????????? ??????????
		Point2Df pos;

		/// ???????????????????? ??????????
		float gam;

	public:

		//#ifndef __CUDACC__
		//		/// MPI-?????????????????? ????????
		//		static MPI_Datatype mpiVortex2Df;
		//#endif
		//
		//		static size_t offsPos;
		//		static size_t offsGam;

				/// ???????????? ??????????????????????
		Vortex2Df() {};

		/// \brief ?????????????????????? ??????????????????????????
		///
		/// \param[in] _r ?????????????????????? ???????????? ???? ????????????-???????????? ?????????????????? ??????????
		/// \param[in] _g ???????????????????? (??????????????????????????) ??????????
		Vortex2Df(const Point2Df& _r, const float _g)
			: pos(_r), gam(_g) { }

		/// ????????????????????
		~Vortex2Df() {};

		/// \brief ?????????????? ?????? ?????????????? ?? ????????????-?????????????? ??????????
		/// \return ???????????? ???? ????????????-???????????? ??????????
		Point2Df& r() { return pos; }

		/// \brief ?????????????? ?????? ?????????????? ?????? ???????????? ?? ????????????-?????????????? ??????????
		/// \return ?????????????????????? ???????????? ???? ????????????-???????????? ??????????
		const Point2Df& r() const { return pos; }

		/// \brief ?????????????? ?????? ?????????????? ?? ???????????????????? ??????????
		/// \return ???????????? ???? ???????????????????? ??????????
		float& g() { return gam; }

		/// \brief ?????????????? ?????? ?????????????? ?????? ???????????? ?? ???????????????????? ??????????
		/// \return ?????????????????????? ???????????? ???? ???????????????????? ??????????
		const float& g() const { return gam; }

		//#ifndef __CUDACC__
		//		/// C?????????????? MPI-?????????????????? ????????
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

	/// ?????????????????????? ???????? ???????????? - ??????????????????, ???????????????? ???? ???? ??????????????????, ?????? ?? ??????????
	typedef Vortex2D Source2D;
	typedef Vortex2Df Source2Df;


}//namespace VMlib

using VMlib::Vortex2D;
using VMlib::Vortex2Df;
#endif
 
