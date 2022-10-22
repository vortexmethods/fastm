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
| File name: operations.cuh                                                   |
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
\brief Вспомогательные операции
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 28 августа 2022 г.
*/

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include <cuda.h>
#include "Params.h"

namespace BHcu
{

#ifndef WIN32
#define __forceinline inline
#endif

	__device__ __forceinline real2 multz(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x * b.x - a.y * b.y;
		res.y = a.x * b.y + a.y * b.x;
		return res;
	}

	__device__ __forceinline real2 multz(real ax, real ay, real2 b)
	{
		real2 res;
		res.x = ax * b.x - ay * b.y;
		res.y = ax * b.y + ay * b.x;
		return res;
	}
	__device__ __forceinline real2 multzA(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x * b.x + a.y * b.y;
		res.y = a.y * b.x - a.x * b.y;
		return res;
	}	
	
	__device__ __forceinline real2 multzA(real2 a, real bx, real by)
	{
		real2 res;
		res.x = a.x * bx + a.y * by;
		res.y = a.y * bx - a.x * by;
		return res;
	}

	__device__ __forceinline real2 operator*(real a, real2 b)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline real2 operator*(real2 b, real a)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline real2 operator/(real2 b, real a)
	{
		real2 res;
		res.x = b.x / a;
		res.y = b.y / a;
		return res;
	}

	__device__ __forceinline real2& operator+=(real2& a, real2 b)
	{
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline real2 operator-(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x - b.x;
		res.y = a.y - b.y;
		return res;
	}

	__device__ __forceinline real2 operator+(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x + b.x;
		res.y = a.y + b.y;
		return res;
	}

}//namespace BHcu

#endif