#ifndef SOURCE_H_
#define SOURCE_H_

#include <vector>
#include <cmath>
#include <string>

#include <omp.h>

#include "numvector.h"

//ВСЕ ВЕКТОРЫ ПОКА ЧТО СДЕЛАНЫ ТРЕХМЕРНЫМИ!!!

#define PI 3.14159265358979323846264338327

inline double dot(const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

//добавил операцию вычисления скалярного произведения
inline double operator* (const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

//Переименовал cross в cross3
inline double cross3(const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	return (x[0] * y[1] - x[1] * y[0]);
}

//Добавил операцию вычисления векторного произведения
inline numvector<double, 3> operator^ (const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	return { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
}

inline double sqr(const double x){ return x*x; }

inline double dist2(const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	return (sqr(x[0] - y[0]) + sqr(x[1] - y[1]));
}

inline numvector<double, 3>& operator*=(numvector<double, 3>& x, const double m)
{
	int dim = x.size();
	for (int i = 0; i < dim; ++i)
		x[i] *= m;
	return x;
}

inline numvector<double, 3> operator+(const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	numvector<double, 3> res(x);
	int dim = x.size();
	for (int i = 0; i < dim; ++i)
		res[i] += y[i];
	return res;
}

inline numvector<double, 3> operator-(const numvector<double, 3>& x, const numvector<double, 3>& y)
{
	numvector<double, 3> res(x);
	int dim = x.size();
	for (int i = 0; i < dim; ++i)
		res[i] -= y[i];
	return res;
}

inline numvector<double, 3> operator*(const double m, const numvector<double, 3>& x)
{
	numvector<double, 3> res(x);
	res *= m;
	return res;
}


void tLayerAver_Matrskos(std::vector<numvector<double, 3>>& CC, std::vector<std::vector<std::vector<double>>>& ME, const std::string& type)
{
	size_t np = CC.size();

	CC.push_back(CC[0]);

	//Panel vectors
	std::vector<numvector<double, 3>> dd;
	dd.resize(np);
	for (size_t i = 0; i < np; ++i)
		dd[i] = { (CC[i + 1][0] - CC[i][0]), (CC[i + 1][1] - CC[i][1]), 0.0 };

	//Lengthes of the panels
	std::vector<double> len;
	len.resize(np);
	for (size_t i = 0; i < np; ++i)
		len[i] = sqrt(dot(dd[i], dd[i]));

	//Tangent unit vectors
	std::vector<numvector<double, 3>> tau(dd);
	for (size_t i = 0; i < np; ++i)
	{
		tau[i] *= 1.0 / len[i];
	};


	if (type == "mtr")
	{
		// Allocate memory for M
		ME.resize(np*2 + 1);
		for (size_t i = 0; i < np*2 + 1; ++i)
		{	
			ME[i].resize(np*2 + 1);  
			for (size_t j = 0; j < np*2 + 1; ++j)
				ME[i][j].resize(1, 0.0);
		}
			
		for (size_t i = 0; i < np; ++i)
		{
			ME[i][np*2][0] = 1.0;     
			ME[np*2][i][0] = len[i]; 
		}
	}

	else
	if (type == "00" || type == "01" || type == "10" || type == "11")	
	{
		ME.resize(np);
		for (size_t i = 0; i < np; ++i)
		{	
			ME[i].resize(np);
			for (size_t j = 0; j < np; ++j)
				ME[i][j].resize(3, 0.0);
		}
	}

		
	//auxillary vectors
	numvector<double, 3> d, d0;

	numvector<double, 3> p1, s1, p2, s2;

    numvector<double, 3> veca, vecaa, vecb, vecbb;

	numvector<double, 3> v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;

	const numvector<double, 3> k ({ 0.0, 0.0, 1.0 });

	//auxillary scalars
	double a[3];
	double qq[3];
	int cft;
	double iznam; 

	int i, j;

	//Собираем все блоки матрицы сразу
	for (i = 0; i < np; ++i)   
	for (j = 0; j < np; ++j)
	{
		//знаменатели у всех блоков одинаковые
		iznam = 1.0 / (sqr(len[i]) * len[j] * PI);

		cft = 1.0;
		if (i != j) //внедиагональные элементы всех блоков
		{
			d = dd[i];
			d0 = dd[j];

			p1 = { (CC[i][0] - CC[j + 1][0]), (CC[i][1] - CC[j + 1][1]), 0.0 };
			s1 = { (CC[i][0] - CC[j][0]), (CC[i][1] - CC[j][1]), 0.0 };
			p2 = { (CC[i + 1][0] - CC[j + 1][0]), (CC[i + 1][1] - CC[j + 1][1]), 0.0 };
			s2 = { (CC[i + 1][0] - CC[j][0]), (CC[i + 1][1] - CC[j][1]), 0.0 };

			if ((j == i + 1) || ((j == 0) && (i == np - 1))) //соседние панели, необходимо перевернуть
			{
				d *= -1.0;
				d0 *= -1.0;
				cft = -1.0;

				p1 = { 0.0, 0.0, 0.0 };
				s1 = { (CC[j][0] - CC[j + 1][0]), (CC[j][1] - CC[j + 1][1]), 0.0 };
				p2 = { (CC[i][0] - CC[j][0]), (CC[i][1] - CC[j][1]), 0.0 };
				s2 = { (CC[i][0] - CC[j + 1][0]), (CC[i][1] - CC[j + 1][1]), 0.0 };

				v5 = { 0.0, 0.0, 0.0 };
				v11 = { 0.0, 0.0, 0.0 };

				a[2] = 0.0;
				qq[2] = 0.0;
			}
			else if ((i == j + 1) || ((i == 0) && (j == np - 1))) //соседние панели, не нужно переворачивать
			{
				v5 = { 0.0, 0.0, 0.0 };
				v11 = { 0.0, 0.0, 0.0 };

				a[2] = 0.0;
				qq[2] = 0.0;
			}
			else //несоседние панели, общий случай
			{
				v5 = (s1*s1)*d0 - (((s1 + s2)*d) / (d*d)) * ((s1*d0)*d + ((s1^d0)^d)) + \
					(d0*d0) / (d*d) * (((s1 + p2)*d)*d + (((s1 + p2)^d)^d));
					
				v11 = (3.0 - 8.0 * (dot(s1, d) * dot(s1, d0)) / (dot(d, d) * dot(d0, d0)) + \
					4.0 * (dot(s1, s1) * dot(d, d0)) / (dot(d, d) * dot(d0, d0)) + 6.0 * dot(d, s1) / dot(d, d) - 6.0 * dot(d0, s1) / dot(d0, d0)) * ((s1*d0)*d + ((s1^d0) ^ d)) + 3.0 * (dot(s1, s1) * p2) - (dot(s1, s1) * s1) - \
					(((d0*d0) / (d*d)) * ((d0*d)*d + ((d0^d) ^ d)));

				a[2] = atan(dot(p2, d) / cross3(p2, d)) - atan(dot(p1, d) / cross3(p1, d));
				
				qq[2] = 0.5 * log(dot(p2, p2) / dot(p1, p1));
			}

			//далее - "универсальные" коэффициенты - для любой конфигурации панелей
			v1 = (d*d) * d0;
			v2 = (s1*d0)*d + ((s1^d0) ^ d);
			v3 = (p1*d0)*d + ((p1^d0) ^ d);
			v4 = (s1*s1)*d0 - (((s1 + s2)*d) / (d*d)) * ((s1*d0)*d + ((s1^d0)^d));
			
			v6 = (d0*d0)*d;
			v7 = (d*d) / (d0*d0) * (((s1 + p2)*d0)*d0 + (((s1 + p2) ^ d0) ^ d0));
			v8 = (2.0* ((p1*d0) / (d0*d0)))*v2 - ((s1*p1)*d + ((s1^p1) ^ d));

			v9 = ((d*d) / (d0*d0))*((d*d0)*d0 + ((d^d0) ^ d0));
			v10 = (3.0 - 8.0 * (dot(s1, d) * dot(s1, d0)) / (dot(d, d) * dot(d0, d0)) + 4.0 * (dot(s1, s1) * dot(d, d0)) / (dot(d, d) * dot(d0, d0)) + 6.0 * dot(d, s1) / dot(d, d) - 6.0 * dot(d0, s1) / dot(d0, d0)) * v2 + 3.0 * (dot(s1, s1) * p2) - (dot(s1, s1) * s1);
			v12 = v6 - v1 - 2.0 * v2;


			a[0] = atan(dot(s2, d0) / cross3(s2, d0)) - atan(dot(p2, d0) / cross3(p2, d0));
			a[1] = atan(dot(s2, d) / cross3(s2, d)) - atan(dot(s1, d) / cross3(s1, d));

			qq[0] = 0.5 * log(dot(s2, s2) / dot(p2, p2));
			qq[1] = 0.5 * log(dot(s2, s2) / dot(s1, s1));

			veca = 0.5 * iznam * ((a[0] * v1 + a[1] * v2 - a[2] * v3) + (k ^ (qq[0] * v1 + qq[1] * v2 - qq[2] * v3)));

			vecaa = 0.25 * iznam * ((a[1] * v4 -  a[2] * v5) + (k ^ (qq[1] * v4 - qq[2] * v5 + v6)));

			vecb = 0.25 * iznam * ((a[0] * v7 + (a[1] - a[2]) * v8) + (k ^ (qq[0] * v7 + (qq[1] - qq[2])*v8 - v1)));

			vecbb = (iznam / 24.0) * ((a[0] * v9 + a[1] * v10 - a[2] * v11) + (k ^ (qq[0] * v9 + qq[1] * v10 - qq[2] * v11 + v12)));

			if (type == "mtr")
			{
				ME[i][j][0] = dot(veca, tau[i]);  //A
				ME[i+np][j][0] = cft*dot(vecaa, tau[i]);  //AA 
				ME[i][j + np][0] = cft * dot(vecb, tau[i]);  //B
				ME[i + np][j + np][0] = dot(vecbb, tau[i]);  //BB  
			}
			else 
			{	
			if (type == "00")
				{
					ME[i][j][0] = veca[1];
					ME[i][j][1] = -veca[0];
					ME[i][j][2] = veca[2];
				}
				else if (type == "01")
				{
					ME[i][j][0] = cft * vecb[1];
					ME[i][j][1] = -cft * vecb[0];
					ME[i][j][2] = cft * vecb[2];
				}
				else if (type == "10")
				{
					ME[i][j][0] = cft * vecaa[1];
					ME[i][j][1] = -cft * vecaa[0];
					ME[i][j][2] = cft * vecaa[2];
				}
				else if (type == "11")
				{
					ME[i][j][0] = vecbb[1];
					ME[i][j][1] = -vecbb[0];
					ME[i][j][2] = vecbb[2];
				}
			}
		}
		else //диагональные элементы всех блоков матрицы
		{
			if (type == "mtr")
			{
				ME[i][j][0] = -0.5;  //A
				ME[i + np][j][0] = 0.0;  //AA
				ME[i][j + np][0] = 0.0;  //B 
				ME[i + np][j + np][0] = -1.0 / 24.0;    //BB 
			}
			if (type == "00" || type == "11")
			{
				ME[i][j][0] = 0.0;
				ME[i][j][1] = 0.0;
				ME[i][j][2] = 0.0;
			}
			if (type == "01")
			{
				numvector<double, 3> d = dd[i];	
				double ln = sqrt(dot(d,d));
				ME[i][j][0] = -0.25/(PI*ln)*d[0];
				ME[i][j][1] = -0.25/(PI*ln)*d[1];
				ME[i][j][2] = -0.25/(PI*ln)*d[2];
			}
			if (type == "10")
			{
				numvector<double, 3> d = dd[i];
				double ln = sqrt(dot(d,d));
				ME[i][j][0] = 0.25/(PI*ln)*d[0];
				ME[i][j][1] = 0.25/(PI*ln)*d[1];
				ME[i][j][2] = 0.25/(PI*ln)*d[2];
			}
		}
	} 
}
#endif