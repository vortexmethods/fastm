/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: defs.h                                                           |
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
\brief Вспомогательные функции
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#pragma once

#include <iostream>
#include "PointsCopy.h"

#ifdef calcOp
	#define ADDOP(n) BH::op += ((int)n)
#else
	#define ADDOP(n) 
#endif

namespace BH
{
	extern long long op;

	//Мировые константы
	static const double PI = 3.1415926535897932384626;
	static const double DPI = 2.0 * 3.1415926535897932384626;
	static const double IPI = 1.0 / PI;
	static const double IDPI = 0.5 / PI;

	/// Длина мортоновского кода для каждой координаты (не более половины длины int в битах)
	static const int codeLength = 14;

	/// 2 в степени длины мортоновского кода (на каждую координату)
	static const int twoPowCodeLength = (1 << codeLength);

	/// \brief Вспомогательная функция, которая определяет, находится ли панель itI на контуре после панели itJ
	///
	/// \param[in] itI константная ссылка на обертку для второй ("правой") панели
	/// \param[in] itJ константная ссылка на обертку для первой ("левой") панели
	inline bool isAfter(const PointsCopy& itI, const PointsCopy& itJ)
	{
		return (itI.panBegin == itJ.panEnd);
	}

	/// \brief Вспомогательная функция вычисления угла между векторами (в диапазоне (-pi...pi]) (со знаком, поворот от первого ко второму)
	///
	/// \param[in] p константная ссылка на первый вектор
	/// \param[in] s константная ссылка на второй вектор 
	inline double Alpha(const Point2D& p, const Point2D& s)
	{
		return atan2(cross3(p, s), p & s);
		ADDOP(5);
	}

	/// \brief  Вспомогательная функция вычисления логарифма отношения норм векторов
	/// 
	/// \param[in] p константная ссылка на первый вектор
	/// \param[in] s константная ссылка на второй вектор 
	inline double Lambda(const Point2D& p, const Point2D& s)
	{
		return 0.5 * log((s & s) / (p & p));
		ADDOP(7);
	}

	/// Вспомогательная функция вычисления влияния панелей, принимает параметры двух панелей
	inline Point2D Lambda(const Point2D& piast, const Point2D& kasiast, double leni, double len1, double s)
	{
		return { (piast[0] + leni * kasiast[0] * s) / len1, (piast[1] + leni * kasiast[1] * s) / len1 };
		ADDOP(6);
	}

	/// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
	inline Point2D Omega(const Point2D& a, const Point2D& b, const Point2D& c)
	{
		return (a & b) * c + (Point2D({ -c[1], c[0] })) * cross3(a, b);
		ADDOP(8);
	}

	/// Вспомогательная функция корректировки capacity вектора (при необходимости - удваивает)
	inline void SizeCheck(std::vector<Point2D>& i00)
	{
		if (i00.capacity() == i00.size())
			i00.reserve(i00.size() * 2);
	}


	/// Умножение комплексных чисел
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
		ADDOP(4);
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[1] * b[0] + a[0] * b[1] });
	}

	/// Возведение в степень комплексных чисел
	inline Point2D powz(const Point2D& z, double n)
	{
		double phi, R;
		ADDOP(10);
		phi = n * atan2(z[1], z[0]);
		R = pow(z.length2(), 0.5*n);
		return Point2D({ R * cos(phi), R * sin(phi) });
	}

	/// Умножение a на комплексно сопряженноe к b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
		ADDOP(4);
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

	/// Шаблонная функция возведения в квадрат
	template <typename T>
	inline T sqr(T x)
	{
		ADDOP(1);
		return x * x;
	}

	/// \brief Шаблонная функция знака числа
	/// Написана оптимизированная версия, которая работает самым быстрым возможным образом, 
	/// т.к. не содержит ни одного оператора условного перехода
	template <typename T>
	inline int sign(T val)
	{
		return (T(0) < val) - (val < T(0));
	}


	//Округление "в потолок" результата деления x на y
	//(годится для любых натуральных, но работает довольно медленно, здесь не нужна)
	//unsigned int ceil(unsigned int x, unsigned int y)
	//{
	//    return x / y + (x % y != 0);
	//}

	/// \brief Округление "в потолок" результата деления на степень двойки, эквивалент ceil(x / (2^p))
	///
	/// \param[in] x делимое
	/// \param[in] p показатель степени двойки в делителе
	inline int ceilpow2(unsigned int x, unsigned int p) 
	{
		return (x >> p) + !!(x & ((1 << p) - 1));
	}

	/// \brief Округление "в потолок" результата деления пополам, эквивалент ceil(x / 2)
	inline int ceilhalf(unsigned int x) 
	{
		return (x >> 1) + (x & 1);
	}

	/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
	///
	/// \param[in] b константная ссылка на вектор или список
	template<typename T>
	inline double norm(const T& b)
	{
		double norm = 0;
#ifndef OLD_OMP
#pragma omp simd reduction(+:norm)
#endif
			for (size_t i = 0; i < b.size(); i++)
				norm += (b[i] * b[i]);
			ADDOP(2*b.size() + 1);
			return sqrt(norm);
	}

	/// Шаблонная функция сложения двух векторов
	template<typename T>
	inline std::vector<T> operator+(const std::vector<T>& x, const std::vector<T>& y)
	{
		std::vector<T> c(x);
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] += y[i];
		return c;
	}


	/// Шаблонная функция прибавления к одному вектору другого
	template<typename T>
	inline std::vector<T>& operator+=(std::vector<T>& x, const std::vector<T>& y)
	{
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			x[i] += y[i];
		return x;
	}

	/// Шаблонная функция вычитания векторов
	template<typename T>
	inline std::vector<T> operator-(const std::vector<T>& x, const std::vector<T>& y)
	{
		std::vector<T> c(x);
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] -= y[i];
		return c;
	}

	/// Шаблонная функция вычитания из одного вектора другого
	template<typename T>
	inline std::vector<T>& operator-=(std::vector<T>& x, const std::vector<T>& y)
	{
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			x[i] -= y[i];
		return x;
	}

	/// Шаблонная функция умножения числа на вектор
	template<typename T>
	inline std::vector<T> operator*(const T lambda, const std::vector<T>& x)
	{
		std::vector<T> c(x);
		c.resize(x.size());
	
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] *= lambda;
		ADDOP(x.size());
		return c;
	}


	/// Шаблонная функция вычисления скалярного произведения двух векторов
	template<typename T>
	inline T operator&(const std::vector<T>& x, const std::vector<T>& y)
	{
		T c = 0;
#ifndef OLD_OMP
#pragma omp simd reduction(+:c)
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c += x[i] * y[i];
		ADDOP(x.size());
		return c;
	}

#ifdef asympScheme
	/// Вспомогательная функцяя для асимптотической схемы
	inline Point2D sFirst(const Point2D& z, const int q, const int p, double len)
	{
		Point2D res = { 0.0, 0.0 };
		Point2D sum, tmp;
		double arg;
		Point2D pw = powz((1.0 / len) * z, -1.0 / q);
		ADDOP(4);

		for (int j = 0; j <= q - 1; ++j)
		{
			arg = DPI * j / q;
			tmp = multz({ cos(arg), sin(arg) }, pw) + Point2D{ 1.0, 0.0 };
			sum = { 0.5 * log(tmp.length2()), atan2(tmp[1], tmp[0]) };
			res += multz({ cos(arg * p), sin(arg * p) }, sum);
			ADDOP(13);
		}
		res = -multz(powz((1.0 / len) * z, (-1.0 * p / q)), res);
		ADDOP(5);

		if (p & 1) //Если p - нечетное
			res *= -1.0;
		ADDOP(2);


		return res;
	}
	inline Point2D Yfrom1(const numvector<Point2D, 2>& cs, const Point2D& idk, const int m, const int j, const int q, const int p)
	{
		Point2D e0, e1, temp, temp2, temp3, log1, mult1, Y;
		double arg, theta;
				temp = multz(cs[m], idk);
				arg = DPI * (j+1) / q;
				e0 = { cos(arg), sin(arg) };
				e1 = { cos(arg * p), sin(arg * p) };
				temp2 = multz(e0, powz(temp, (-1.0 / q)));
				temp3 = { temp2[0] + 1, temp2[1] };
				log1 = { 0.5 * log(temp3.length2()), atan2(temp3[1], temp3[0]) };
				mult1 = multz(powz(temp, (double)(1.0 -(1.0 * p / q))), log1);
				Y = multz(mult1, e1);
		return multz(mult1, e1);
	}

	inline numvector<Point2D, 2> StoConstFrom1_new(const numvector<Point2D, 2>& pan1arb, const numvector<Point2D, 2>& paniarb, const int t, const params& prm)
	{
		auto q = prm.q;
		auto p = prm.p;
		Point2D e, temp, temp1, temp2, temp3, temp4, H0, H1, log1, mult1, mult2, sum, sum1, sum2, sum3, H0_math, H1_math, h;
		double arg, theta;
		sum = { 0.0, 0.0 };
		sum3 = { 0.0, 0.0 };
		Point2D dk = (pan1arb[1] - pan1arb[0]);
		Point2D di = (paniarb[1] - paniarb[0]);
		Point2D idk = (1.0 / (dk.length2())) * Point2D({ dk[0], -dk[1] });
		Point2D idi = (1.0 / (di.length2())) * Point2D({ di[0], -di[1] });
		numvector<Point2D, 2> cs = { (paniarb[1] - pan1arb[0]), (paniarb[0] - pan1arb[0]) };
		numvector<Point2D, 2> cp = { (paniarb[1] - pan1arb[1]), (paniarb[0] - pan1arb[1]) }; 
		Point2D icp1 = (1.0 / (cp[1].length2())) * Point2D({cp[1][0], -cp[1][1]});
		double len1 = sqrt((pan1arb[1] - pan1arb[0]) & (pan1arb[1] - pan1arb[0]));
		double leni = sqrt((paniarb[1] - paniarb[0]) & (paniarb[1] - paniarb[0]));
		Point2D kas1 = { (pan1arb[1] - pan1arb[0])[0] / len1, (pan1arb[1] - pan1arb[0])[1] / len1 };

		//Поиск H0
			theta = atan2(kas1[1], kas1[0]);
			e = { cos(theta), -sin(theta) };
			mult1 = ( leni/(DPI * (1.0 - (1.0 * p[t] / q[t]))) ) * e;
			mult2 = multz(dk, idi);
			temp1 = multz(cp[0], icp1);
			log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };
			for (int j = 0; j < q[t]; j++) {
				sum += Yfrom1(cs, idk, 0, j, q[t], p[t]) - Yfrom1(cs, idk, 1, j, q[t], p[t]);
			}
			if (p[t] & 1)
				sum *= -1.0;
			temp2 = multz(mult1, mult2);
			H0_math = multz(temp2, log1 - sum);
			H0 = { H0_math[0], -H0_math[1] };

		//Поиск H1
			theta = atan2(kas1[1], kas1[0]);
			e = { cos(theta), -sin(theta) };
			mult1 = ( leni / (DPI * (2.0 - (1.0 * p[t] / q[t]))) ) * e;
			mult2 = powz(multz(dk, idi), 2);
			temp1 = multz(cp[0], icp1);
			log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };

			temp = (2.0 - (1.0 * p[t] / q[t])) / (1.0 - (1.0 * p[t] / q[t])) * multz(cs[1], idk);
			h = { 1 - temp[0], -temp[1] };

			sum1 = multz(h, log1);
			sum2 = multz(di, idk);
			temp2 = h + multz(cp[0], idk);
			temp3 = h + multz(cp[1], idk);
			for (int j = 0; j < q[t]; j++) {
				sum3 += multz(temp2, Yfrom1(cs, idk, 0, j, q[t], p[t])) - multz(temp3, Yfrom1(cs, idk, 1, j, q[t], p[t]));
			}
			if (p[t] & 1)
				sum3 *= -1.0;
			temp4 = multz(mult1, mult2);
			H1_math = multz(temp4, sum1 + sum2 - sum3); 
			//H1_math = multz(temp4, sum1 + sum2 - sum3) - 0.5 * H0;
			H1 = { H1_math[0], -H1_math[1] };

		return { H0, H1 };
	}

	inline numvector<Point2D, 2> StoConstFrom1_new1(const numvector<Point2D, 2>& pan1arb, const numvector<Point2D, 2>& paniarb, const int t, const params& prm)
	{
		auto q = prm.q;
		auto p = prm.p;
		Point2D e, temp, temp1, temp2, temp3, temp4, H0, H1, log1, mult1, mult2, sum, sum1, sum2, sum3, H0_math, H1_math, h, sum4;
		double arg, theta;
		sum = { 0.0, 0.0 };
		sum1 = { 0.0, 0.0 };
		sum3 = { 0.0, 0.0 };
		sum4 = { 0.0, 0.0 };
		Point2D dk = (pan1arb[1] - pan1arb[0]);
		Point2D di = (paniarb[1] - paniarb[0]);
		Point2D idk = (1.0 / (dk.length2())) * Point2D({ dk[0], -dk[1] });
		Point2D idi = (1.0 / (di.length2())) * Point2D({ di[0], -di[1] });
		numvector<Point2D, 2> cs = { (paniarb[1] - pan1arb[0]), (paniarb[0] - pan1arb[0]) };
		numvector<Point2D, 2> cp = { (paniarb[1] - pan1arb[1]), (paniarb[0] - pan1arb[1]) };
		Point2D icp1 = (1.0 / (cp[1].length2())) * Point2D({ cp[1][0], -cp[1][1] });
		double len1 = sqrt((pan1arb[1] - pan1arb[0]) & (pan1arb[1] - pan1arb[0]));
		double leni = sqrt((paniarb[1] - paniarb[0]) & (paniarb[1] - paniarb[0]));
		Point2D kas1 = { (pan1arb[1] - pan1arb[0])[0] / len1, (pan1arb[1] - pan1arb[0])[1] / len1 };

		//Поиск H0
		theta = atan2(kas1[1], kas1[0]);
		e = { cos(theta), -sin(theta) };
		mult1 = (leni / (DPI * (1.0 - (1.0 * p[t] / q[t])))) * e;
		mult2 = multz(dk, idi);
		temp1 = (1.0 / q[t]) * multz(di, idk);
		log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };
		for (int j = 0; j < q[t]; j++) {
			sum += Yfrom1(cs, idk, 0, j, q[t], p[t]);
			if((j+1) !=  q[t] / 2)
			sum4 += Yfrom1(cs, idk, 1, j, q[t], p[t]);
		}
		sum1 = sum - sum4;
		if (p[t] & 1)
			sum1 *= -1.0;
		temp2 = multz(mult1, mult2);
		H0_math = multz(temp2, log1 - sum1);
		H0 = { H0_math[0], -H0_math[1] };

		//Поиск H1
		theta = atan2(kas1[1], kas1[0]);
		e = { cos(theta), -sin(theta) };
		mult1 = (leni / (DPI * (2.0 - (1.0 * p[t] / q[t])))) * e;
		mult2 = powz(multz(dk, idi), 2);
		temp1 = (1.0 / q[t]) * multz(di, idk);
		log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };

		temp = (2.0 - (1.0 * p[t] / q[t])) / (1.0 - (1.0 * p[t] / q[t])) * multz(cs[1], idk);
		h = { 1 - temp[0], -temp[1] };

		sum1 = multz(h, log1);
		sum2 = multz(di, idk);
		temp2 = h + multz(cp[0], idk);
		temp3 = h;
		for (int j = 0; j < q[t]; j++) {
			sum3 += multz(temp2, Yfrom1(cs, idk, 0, j, q[t], p[t]));
			if ((j+1) != q[t]/2)
				sum3 -= multz(temp3, Yfrom1(cs, idk, 1, j, q[t], p[t]));
		}
		if (p[t] & 1)
			sum3 *= -1.0;
		temp4 = multz(mult1, mult2);
		H1_math = multz(temp4, sum1 + sum2 - sum3);
		H1 = { H1_math[0], -H1_math[1] };
		//H1 = Point2D({ H1_math[0], -H1_math[1] }) - 0.5 * H0;
		return { H0, H1 };
	}

	inline numvector<Point2D, 2> StoConstFrom1_new2(const numvector<Point2D, 2>& pan1arb, const numvector<Point2D, 2>& paniarb, const int t, const params& prm)
	{
		auto q = prm.q;
		auto p = prm.p;
		Point2D e, temp, temp1, temp2, temp3, temp4, H0, H1, log1, mult1, mult2, sum, sum1, sum2, sum3, H0_math, H1_math, h;
		double arg, theta;
		sum = { 0.0, 0.0 };
		sum3 = { 0.0, 0.0 };
		Point2D dk = (pan1arb[1] - pan1arb[0]);
		Point2D di = (paniarb[1] - paniarb[0]);
		Point2D idk = (1.0 / (dk.length2())) * Point2D({ dk[0], -dk[1] });
		Point2D idi = (1.0 / (di.length2())) * Point2D({ di[0], -di[1] });
		numvector<Point2D, 2> cs = { (paniarb[1] - pan1arb[0]), (paniarb[0] - pan1arb[0]) };
		numvector<Point2D, 2> cp = { (paniarb[1] - pan1arb[1]), (paniarb[0] - pan1arb[1]) };
		Point2D icp1 = (1.0 / (cp[1].length2())) * Point2D({ cp[1][0], -cp[1][1] });
		double len1 = sqrt((pan1arb[1] - pan1arb[0]) & (pan1arb[1] - pan1arb[0]));
		double leni = sqrt((paniarb[1] - paniarb[0]) & (paniarb[1] - paniarb[0]));
		Point2D kas1 = { (pan1arb[1] - pan1arb[0])[0] / len1, (pan1arb[1] - pan1arb[0])[1] / len1 };

		//Поиск H0
		theta = atan2(kas1[1], kas1[0]);
		e = { cos(theta), -sin(theta) };
		mult1 = (leni / (DPI * (1.0 - (1.0 * p[t] / q[t])))) * e;
		mult2 = multz(dk, idi);
		temp1 = multz(cp[0], icp1);
		log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };
		for (int j = 0; j < q[t]; j++) {
			sum += Yfrom1(cs, idk, 1, j, q[t], p[t]);
		}
		if (p[t] & 1)
			sum *= -1.0;
		temp2 = multz(mult1, mult2);
		H0_math = multz(temp2, log1 + sum);
		H0 = { H0_math[0], -H0_math[1] };

		//Поиск H1
		theta = atan2(kas1[1], kas1[0]);
		e = { cos(theta), -sin(theta) };
		mult1 = (leni / (DPI * (2.0 - (1.0 * p[t] / q[t])))) * e;
		mult2 = powz(multz(dk, idi), 2);
		temp1 = multz(cp[0], icp1);
		log1 = { 0.5 * log(temp1.length2()), atan2(temp1[1], temp1[0]) };

		temp = (2.0 - (1.0 * p[t] / q[t])) / (1.0 - (1.0 * p[t] / q[t])) * multz(cs[1], idk);
		h = { 1 - temp[0], -temp[1] };

		sum1 = multz(h, log1);
		sum2 = multz(di, idk);
		temp3 = h + multz(cp[1], idk);
		for (int j = 0; j < q[t]; j++) {
			sum3 += multz(temp3, Yfrom1(cs, idk, 1, j, q[t], p[t]));
		}
		if (p[t] & 1)
			sum3 *= -1.0;
		temp4 = multz(mult1, mult2);
		H1_math = multz(temp4, sum1 + sum2 + sum3);
		H1 = { H1_math[0], -H1_math[1] };
		//H1 = Point2D({ H1_math[0], -H1_math[1] }) - 0.5 * H0;
		return { H0, H1 };
	}
#endif

}//namespace BH
