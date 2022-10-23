/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.2    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/10/22     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
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
\version 1.2
\date 22 октября 2022 г.
*/

#pragma once

#include <iostream>
#include "PointsCopy.h"

namespace BH
{
	// Вспомогательная функция, которая определяет, находится ли панель itI после панели itJ
	inline bool isAfter(const PointsCopy& itI, const PointsCopy& itJ)
	{
		return (itI.panBegin == itJ.panEnd);
	}

	// Вспомогательная функция вычисления угла между векторами
	inline double Alpha(const Point2D& p, const Point2D& s)
	{
		return atan2(cross3(p, s), p & s);
	}

	// Вспомогательная функция вычисления логарифма отношения норм векторов
	inline double Lambda(const Point2D& p, const Point2D& s)
	{
		return 0.5 * log((s & s) / (p & p));
	}

	// Вспомогательная функция вычисления 
	inline Point2D Lambda(const Point2D& piast, const Point2D& kasiast, double leni, double len1, double s)
	{
		return { (piast[0] + leni * kasiast[0] * s) / len1, (piast[1] + leni * kasiast[1] * s) / len1 };
	}

	// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
	inline Point2D Omega(const Point2D& a, const Point2D& b, const Point2D& c)
	{
		return (a & b) * c + (Point2D({ -c[1], c[0] })) * cross3(a, b);
	}

	// Вспомогательная функция корректировки capacity вектора
	inline void SizeCheck(std::vector<Point2D>& i00)
	{
		if (i00.capacity() == i00.size())
			i00.reserve(i00.size() * 2);
	}


	//умножение комплексных чисел
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[1] * b[0] + a[0] * b[1] });
	}

	//возведение в степень комплексных чисел
	inline Point2D powz(const Point2D& z, double n)
	{
		double phi, R;
#ifdef calcOp
		op += 4;
#endif
		phi = n * atan2(z[1], z[0]);
		R = pow(z.length2(), 0.5*n);
		return Point2D({ R * cos(phi), R * sin(phi) });
	}

	// умножение a на комплексно сопряженноe к b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

	//возведение в квадрат
	template <typename T>
	inline T sqr(T x)
	{
		return x * x;
	}

	//Знак числа, написана оптимизированная версия, 
	//которая работает самым быстрым возможным образом, т.к. не содержит ни одного оператора условного перехода
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

	//Округление "в потолок" результата деления на степень двойки
	inline int ceilpow2(unsigned int x, unsigned int p) //  =ceil(x / 2^p)
	{
		return (x >> p) + !!(x & ((1 << p) - 1));
	}

	//Округление "в потолок" результата деления пополам
	inline int ceilhalf(unsigned int x) //  =ceil(x / 2), т.е. предыдущая функция при p=1
	{
		return (x >> 1) + (x & 1);
	}

	// Нормы
	template< typename T>
	inline double norm(const unsigned int DIM, const  T& b, const char flag)
	{
		double  norm = 0;
		if (flag == 'k')  //кубическая норма
		{
			for (unsigned int i = 0; i < DIM; i++)
				if (norm < fabs(b[i]))  norm = fabs(b[i]);
			return norm;
		}

		if (flag == '1') //октаэдрическая норма
		{
			for (unsigned int j = 0; j < DIM; j++)
				norm += fabs(b[j]);
			return norm;
		}

		if (flag == '2')  //eвклидова
		{
			for (unsigned int i = 0; i < DIM; i++)
				norm += (b[i] * b[i]);
			return sqrt(norm);
		}
		return 0;
	}

	// return x + y
	inline std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& y)
	{
		std::vector<double> c(x);
		for (int i = 0; i < (int)x.size(); ++i)
			c[i] += y[i];
		return c;
	}

	// x += y
	inline void operator+=(std::vector<double>& x, const std::vector<double>& y)
	{
		for (size_t i = 0; i < x.size(); ++i)
			x[i] += y[i];
	}

	// return x - y
	inline std::vector<double> operator-(const std::vector<double>& x, const std::vector<double>& y)
	{
		std::vector<double> c(x);
		for (int i = 0; i < (int)x.size(); ++i)
			c[i] -= y[i];
		return c;
	}

	// x -= y
	inline void operator-=(std::vector<double>& x, const std::vector<double>& y)
	{
		for (size_t i = 0; i < x.size(); ++i)
			x[i] -= y[i];
	}

	// return lambda * x 
	inline std::vector<double> operator*(const double lambda, const std::vector<double>& x)
	{
		std::vector<double> c;
		c.resize(x.size());
		for (size_t i = 0; i < x.size(); ++i)
			c[i] = lambda * x[i];
		return c;
	}

	// return (x, y) 
	inline double operator&(const std::vector<double>& x, const std::vector<double>& y)
	{
		double c = 0.0;
		for (size_t i = 0; i < x.size(); ++i)
			c += x[i] * y[i];
		return c;
	}

	inline double norm(const std::vector<double>& x1, const std::vector<double>& x2)
	{
		std::vector <double> kk(x1.size());
		for (size_t i = 0; i < x1.size(); ++i)
			kk[i] = std::fabs(x1[i] - x2[i]);

		double res;
		res = kk[0];
		for (size_t i = 1; i < x1.size(); ++i)
		{
			if (res < kk[i])
				res = kk[i];
		}
		return res;
	}

	inline Point2D sFirst(const Point2D& z, const int q, const int p, double len)
	{
		Point2D res = { 0.0, 0.0 };
		Point2D sum, tmp;
		double arg;
		Point2D pw = powz((1.0 / len) * z, -1.0 / q);
		for (int j = 0; j <= q - 1; ++j)
		{
			arg = 2.0 * PI * j / q;
			tmp = multz({ cos(arg), sin(arg) }, pw) + Point2D{ 1.0, 0.0 };
			sum = { 0.5 * log(tmp.length2()), atan2(tmp[1], tmp[0]) };
			res += multz({ cos(arg * p), sin(arg * p) }, sum);
		}
		res = -multz(powz((1.0 / len) * z, (-1.0 * p / q)), res);

		if (p & 1) //Если p - нечетное
			res *= -1.0;

		return res;
	}

	// Вспомогательная функция вычисления 
	inline Point2D F(const Point2D& z, const int q, const int p)
	{
		Point2D res = { 0.0, 0.0 };
		Point2D e, sum, sum1, sum2;
		double arg;
		Point2D pw = powz(z, -1.0 / q);
		for (int j = 0; j <= q - 1; ++j)
		{
			arg = 2.0 * PI * j / q;
			e = { cos(arg * p), sin(arg * p) };
			sum1 = multz({ cos(arg), sin(arg) }, pw);
			sum2 = { sum1[0] + 1.0, sum1[1] };
			sum = { 0.5 * log(sum2.length2()), atan2(sum2[1], sum2[0]) };
			res += multz(e, sum);
		}

		if (p & 1) //Если p - нечетное
			res *= -1.0;

		return res;
	}

	// Вспомогательная функция вычисления 
	inline Point2D Ftilde(const Point2D& z, const int q, const int p)
	{
		Point2D res = { 0.0, 0.0 };
		Point2D e, sum, sum1, sum2;
		double arg;
		Point2D pw = powz(z, -1.0 / q);
		for (int j = 0; j <= q / 2 - 1; ++j)
		{
			arg = 2.0 * PI * j / q;
			e = { cos(arg * p), sin(arg * p) };
			sum1 = multz({ cos(arg), sin(arg) }, pw);
			sum2 = { sum1[0] + 1, sum1[1] };
			sum = { 0.5 * log(sum2.length2()), atan2(sum2[1], sum2[0]) };
			res += multz(e, sum);
		}
		for (int j = q / 2 + 1; j <= q - 1; ++j)
		{
			arg = 2.0 * PI * j / q;
			e = { cos(arg * p), sin(arg * p) };
			sum1 = multz({ cos(arg), sin(arg) }, pw);
			sum2 = { sum1[0] + 1, sum1[1] };
			sum = { 0.5 * log(sum2.length2()), atan2(sum2[1], sum2[0]) };
			res += multz(e, sum);
		}

		if (p & 1)
			res *= -1.0;

		return res;
	}

	// Вспомогательная функция вычисления 
	inline Point2D Neopr0(const Point2D& piast, const Point2D& kasiast, double leni, double len1, double s, const int q, const int p)
	{
		Point2D lam = Lambda(piast, kasiast, leni, len1, s);
		Point2D sum1 = { 0.5 * log((lam[0] - 1) * (lam[0] - 1) + (lam[1] * lam[1])), atan2(lam[1], lam[0] - 1) };
		Point2D sum2 = multz(powz(lam, (double)(q - p) / q), F(lam, q, p));
		Point2D icmplx = (1.0 / (kasiast.length2())) * Point2D({ kasiast[0], -kasiast[1] });

		return (q * len1 / (2.0 * PI * (q - p))) * multz(icmplx, sum1 - sum2);
	}

	// Вспомогательная функция вычисления 
	inline Point2D NeoprLow0(const Point2D& kasiast, double len1, double s, const int q, const int p)
	{
		Point2D sum1 = { log(q), 0.0 };
		Point2D sum2 = Ftilde({ 1.0, 0.0 }, q, p);
		Point2D icmplx = (1.0 / (kasiast.length2())) * Point2D({ kasiast[0], -kasiast[1] });

		return (q * len1 / (2.0 * PI * (q - p))) * multz(icmplx, sum1 - sum2);
	}

	// Вспомогательная функция вычисления 
	inline Point2D H(const Point2D& piast, double len1, const int q, const int p)
	{
		return { 1.0 - (2.0 * q - p) * piast[0] / ((q - p)  * len1),  -((2.0 * q - p) * piast[1] / ((q - p) * len1)) };
	}

	// Вспомогательная функция вычисления 
	inline Point2D Neopr1(const Point2D& piast, const Point2D& kasiast, double leni, double len1, double s, const int q, const int p)
	{
		Point2D h = H(piast, len1, q, p);
		Point2D lam = Lambda(piast, kasiast, leni, len1, s);
		Point2D lam1 = { lam[0] - 1, lam[1] };
		Point2D sum = multz(powz(lam, (double)(q - p) / q), F(lam, q, p));
		Point2D sum1 = { 0.5 * log(lam1.length2()), atan2(lam1[1], lam1[0]) };
		Point2D sum2 = { sum[0] - 1, sum[1] };

		Point2D kasiast2 = multz(kasiast, kasiast);
		Point2D icmplx = (1.0 / (kasiast2.length2())) * Point2D({ kasiast2[0], -kasiast2[1] });

		Point2D sum3 = { multz(h, sum1)[0] - (h[0] - 1) - multz((h + lam1), sum2)[0], multz(h, sum1)[1] - (h[1]) - multz((h + lam1), sum2)[1] };

		return 1.0 * q * len1 * len1 / (2.0 * PI * (2.0 * q - p) * leni) * multz(icmplx, sum3);
	}

	// Вспомогательная функция вычисления 
	inline Point2D NeoprLow1(const Point2D& kasiast, double leni, double len1, const int q, const int p)
	{
		Point2D sum1 = { log(q), 0.0 };
		Point2D sum2 = Ftilde({ 1.0, 0.0 }, q, p);
		Point2D kasiast2 = multz(kasiast, kasiast);
		Point2D icmplx = (1.0 / (kasiast2.length2())) * Point2D({ kasiast2[0], -kasiast2[1] });
		Point2D sum = { 1.0 - q / (q - p) * (sum1[0] - sum2[0]), q / (q - p) * (sum1[1] - sum2[1]) };

		return 1.0 * q * len1 * len1 / (2.0 * PI * (2.0 * q - p) * leni) * multz(icmplx, sum);
	}

#ifdef asympScheme
	// Вспомогательная функция вычисления вектора
	inline numvector<Point2D, 2> StoConstFrom1(const numvector<Point2D, 2>& pan1arb, const numvector<Point2D, 2>& paniarb, const int t)
	{
		double len1 = sqrt((pan1arb[1] - pan1arb[0]) & (pan1arb[1] - pan1arb[0]));
		double leni = sqrt((paniarb[1] - paniarb[0]) & (paniarb[1] - paniarb[0]));
		Point2D kas1 = { (pan1arb[1] - pan1arb[0])[0] / len1, (pan1arb[1] - pan1arb[0])[1] / len1 };
		Point2D kasi = { (paniarb[1] - paniarb[0])[0] / leni, (paniarb[1] - paniarb[0])[1] / leni };
		Point2D nrm1 = { kas1[1], -kas1[0] };
		Point2D nrmi = { kasi[1], -kasi[0] };
		Point2D axx = kas1;
		Point2D ayy = -nrm1;
		Point2D nrmiast = { nrmi & axx, nrmi & ayy };
		Point2D kasiast = { kasi & axx, kasi & ayy };
		Point2D piast = { (paniarb[0] - pan1arb[0]) & axx, (paniarb[0] - pan1arb[0]) & ayy };
		Point2D piendast = { (paniarb[1] - pan1arb[0]) & axx, (paniarb[1] - pan1arb[0]) & ayy };
		Point2D Delta0, Delta1, Delta;
		Point2D h = H(piast, len1, q[t], p[t]);
		Point2D lam0 = Lambda(piast, kasiast, leni, len1, 0);
		Point2D lam1 = Lambda(piast, kasiast, leni, len1, 1);
		Point2D lower0, lower1, upper0, upper1;
		Point2D S00, S0, S11, S1;

		Point2D icmplx = (1.0 / (kasiast.length2())) * Point2D({ kasiast[0], -kasiast[1] });
		Point2D kasiast2 = multz(kasiast, kasiast);
		Point2D icmplx2 = (1.0 / (kasiast2.length2())) * Point2D({ kasiast2[0], -kasiast2[1] });

		if ((piast[1] * piendast[1] > 0) || (piast[0] > 0))
			Delta0 = { 0.0, 0.0 };
		else
		{
			if (piast[1] > 0)
				Delta0 = { -q[t] * len1 / (q[t] - p[t]) * icmplx[0], q[t] * len1 / (q[t] - p[t]) * icmplx[1] };
			else
				Delta0 = { q[t] * len1 / (q[t] - p[t]) * icmplx[0], -q[t] * len1 / (q[t] - p[t]) * icmplx[1] };
		}

		if ((piast[1] * piendast[1] > 0) || (piast[0] > 0))
			Delta1 = { 0.0, 0.0 };
		else
		{
			if (piast[1] > 0) {
				Delta1 = q[t] * len1 * len1 / ((2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h);
			}
			else
				Delta1 = -q[t] * len1 * len1 / ((2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h);;
		}

		if (sqrt((lam0[0] - 1) * (lam0[0] - 1) + lam0[1] * lam0[1]) < 1e-6)
		{
			Delta0 = { 0.0, 0.0 };
			lower0 = NeoprLow0(kasiast, len1, 0, q[t], p[t]);
		}
		else
			lower0 = Neopr0(piast, kasiast, leni, len1, 0, q[t], p[t]);

		if (sqrt((lam0[0] - 1) * (lam0[0] - 1) + lam0[1] * lam0[1]) <= 1e-6)
		{
			Delta1 = { 0.0, 0.0 };
			lower1 = NeoprLow1(kasiast, leni, len1, q[t], p[t]);
		}
		else
			lower1 = Neopr1(piast, kasiast, leni, len1, 0, q[t], p[t]);


		if (sqrt(lam1.length2()) < 1e-6)
		{
			if (piast[1] > 0)
				Delta0 = { -q[t] * len1 / (2.0 * (q[t] - p[t])) * icmplx[1], q[t] * len1 / (2.0 * (q[t] - p[t])) * icmplx[0] };
			else
				Delta0 = { q[t] * len1 / (2.0 * (q[t] - p[t])) * icmplx[1], -q[t] * len1 / (2.0 * (q[t] - p[t])) * icmplx[0] };
			upper0 = { 0.0, 0.0 };
		}
		else
			upper0 = Neopr0(piast, kasiast, leni, len1, 1, q[t], p[t]);

		if (sqrt(lam1.length2()) < 1e-6)
		{
			if (piast[1] > 0)
				Delta1 = { -q[t] * len1 * len1 / (2.0 * (2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h)[1], q[t] * len1 * len1 / (2.0 * (2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h)[0] };
			else
				Delta1 = { q[t] * len1 * len1 / (2.0 * (2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h)[1], -q[t] * len1 * len1 / (2.0 * (2.0 * q[t] - p[t]) * leni) * multz(icmplx2, h)[0] };
			upper1 = { 0.0, 0.0 };
		}
		else
			upper1 = Neopr1(piast, kasiast, leni, len1, 1, q[t], p[t]);

		S00 = { upper0[0] - lower0[0] + Delta0[0], -upper0[1] + lower0[1] - Delta0[1] };
		S0 = { axx * S00[0] + ayy * S00[1] };

		S11 = { upper1[0] - lower1[0] + Delta1[0], -upper1[1] + lower1[1] - Delta1[1] };
		S1 = { axx * S11[0] + ayy * S11[1] };

		return { S0, S1 };
	}
#endif

}//namespace BH
