/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Methods.h                                                        |
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
\brief Процедуры итерационных решателей
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/

#include <iostream>
#include "BarnesHut.h"

namespace BH
{

	const double epsNorm = 1e-10;

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
		for (int i = 0; i < x.size(); ++i)
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
		for (int i = 0; i < x.size(); ++i)
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

	/*
	void BiCGStab(BarnesHut& BH, std::vector<double>& X, double R, std::vector<double>& rhs, std::vector<double>& len, std::vector<Point2D>& tau, int n, int type)
	{
		std::vector<Point2D> AX(n);

		std::vector<double> r(n + 1), rconst(n + 1), p(n + 1), s(n + 1), Ap(n + 1), As(n + 1);

		double timing2 = 0.0, timing3 = 0.0;

		BH.InfluenceComputation(AX, type, timing2, timing3);
		for (int i = 0; i < n; ++i)
		{
			r[i] = rhs[i] - (AX[i] & tau[i]);
			rconst[i] = r[i];
			p[i] = r[i];
		}
		r[n] = len & X;
		rconst[n] = r[n];
		p[n] = r[n];

		double alpha, omega, den, beta;

		std::vector <double> id2len(n);
		for (size_t i = 0; i < len.size(); ++i)
			id2len[i] = 1.0 / (2.0 * len[i]);

		int it = 0;
		do
		{
			BH.IterativeInfluenceComputation(AX, p, type);

			for (int i = 0; i < n; ++i)
				Ap[i] = (AX[i] & tau[i]) + p[n] - p[i] * id2len[i];
			Ap[n] = len & p;

			alpha = (rconst & r) / (Ap & rconst);

			s = r - (alpha * Ap);

			BH.IterativeInfluenceComputation(AX, s, type);


			for (int i = 0; i < n; ++i)
				As[i] = (AX[i] & tau[i]) + s[n] - s[i] * id2len[i];
			As[n] = len & s;

			omega = (As & s) / (As & As);

			X += (alpha * p) + (omega * s);
			R += alpha * p[n] + omega * s[n];

			den = omega * (r & rconst);

			r = s - (omega * As);
			beta = alpha * (r & rconst) / den;

			p -= (omega * Ap);
			p = r + (beta * p);

			it++;
		} while (norm(n, r, '2') > epsNorm);

		std::cout << "BiCGStab: " << it << std::endl;

	}
	*/

	bool IterRot(const std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, double& c, double& s, int m, int n)
	{
		bool fl;
		double hii;

		if (m == 1)
			hii = H[0][0];
		else
			hii = c * H[m - 1][m - 1] - s * H[m - 2][m - 1];

		c = hii / sqrt(hii * hii +
			H[m][m - 1] * H[m][m - 1]);
		s = H[m][m - 1] / sqrt(hii * hii +
			H[m][m - 1] * H[m][m - 1]);

		gs = -s * gs;

		if ((fabs(gs) / norm(n, rhs, '2')) < epsNorm) fl = true;
		else fl = false;
		return fl;
	}

	void SolCircleRun(std::vector<double>& AX, const std::vector<double>& len, const std::vector<Point2D>& tau, const std::vector<double>& rhs, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();

		std::vector<double> alpha(n), beta(n), delta(n), phi(n), xi(n);
		double yn, a;

		alpha[1] = (tau[0] & pnt[0].c1) * len[0] * 24.0;
		beta[1] = -(tau[n - 1] & pnt[0].a1) * len[0] * 24.0; ///почему тут минус?
		delta[1] = -24.0 * len[0] * rhs[n + 0];
		double zn;

		for (int i = 1; i < n - 1; ++i)
		{
			a = (tau[i] & pnt[i].a1);
			zn = a * alpha[i] - 1.0 / (24.0 * len[i]);
			alpha[i + 1] = -(tau[i] & pnt[i].c1) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			delta[i + 1] = (rhs[n + i] - a * delta[i]) / zn;
		}

		a = tau[n - 2] & pnt[n - 2].a1;
		zn = alpha[n - 2] * a - 1.0 / (24.0 * len[n - 2]);
		phi[n - 1] = -((tau[n - 2] & pnt[n - 2].c1) + beta[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n + n - 2] - delta[n - 2] * a) / zn;

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
		}

		double e = (tau[n - 1] & pnt[n - 1].c1);
		a = (tau[n - 1] & pnt[n - 1].a1);

		yn = (rhs[n + n - 1] - e * xi[1] - a * xi[n - 1]) / (e * phi[1] + a * phi[n - 1] - 0.5 / len[n - 1]);

		AX[n + n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
			AX[n + i] = phi[i + 1] * yn + xi[i + 1];

	}

	void SolM(std::vector<double>& AX, const std::vector<double>& len, const std::vector<Point2D>& tau, const std::vector<double>& rhs, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();

		std::vector<double> alpha(n), beta(n), gamma(n), delta(n), phi(n), xi(n), psi(n);
		double yn, ynn;
		double lam1, lam2, mu1, mu2, xi1, xi2;
		double a;

		alpha[1] = (tau[0] & pnt[0].c) * len[0] * 2.0;
		beta[1] = -(tau[n - 1] & pnt[0].a) * len[0] * 2.0; ///почему тут минус?
		gamma[1] = 2.0 * len[0];
		delta[1] = -2.0 * len[0] * rhs[0];
		double zn;

		for (int i = 1; i < n - 1; ++i)
		{
			a = (tau[i] & pnt[i].a);
			zn = a * alpha[i] - 0.5 / len[i];
			alpha[i + 1] = -(tau[i] & pnt[i].c) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			gamma[i + 1] = -(1.0 + a * gamma[i]) / zn;
			delta[i + 1] = (rhs[i] - a * delta[i]) / zn;
		}

		a = tau[n - 2] & pnt[n - 2].a;
		zn = alpha[n - 2] * a - 0.5 / len[n - 2];
		phi[n - 1] = -((tau[n - 2] & pnt[n - 2].c) + beta[n - 2] * a) / zn;
		psi[n - 1] = -(1.0 + gamma[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n - 2] - delta[n - 2] * a) / zn;

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			psi[i] = alpha[i] * psi[i + 1] + gamma[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
		}
		double e = (tau[n - 1] & pnt[n - 1].c);
		a = (tau[n - 1] & pnt[n - 1].a);
		lam1 = e * phi[1] + a * phi[n - 1] - 0.5 / len[n - 1];
		mu1 = e * psi[1] + a * psi[n - 1] + 1.0;
		xi1 = rhs[n - 1] - e * xi[1] - a * xi[n - 1];
		lam2 = mu2 = xi2 = 0.0;
		for (int j = 0; j < n - 1; ++j)
		{
			lam2 += phi[j + 1];
			mu2 += psi[j + 1];
			xi2 -= xi[j + 1];
		}
		lam2 += 1.0;
#ifndef linScheme
		xi2 += rhs[n];
#else
		xi2 += rhs[2 * n];
#endif
		zn = lam1 * mu2 - lam2 * mu1;
		yn = (xi1 * mu2 - xi2 * mu1) / zn;
		ynn = -(xi1 * lam2 - xi2 * lam1) / zn;

#ifndef linScheme
		AX[n] = ynn;
#else
		AX[2 * n] = ynn;
#endif

		AX[n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
			AX[i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];

		//Циклическая прогонка для линейной схемы
#ifdef linScheme 
		SolCircleRun(AX, len, tau, rhs, pnt);
#endif 
	}


	/*
	void GMRES(BarnesHut& BH, std::vector <double>& X, double R, std::vector<double>& rhs, std::vector<double>& len, std::vector<Point2D>& tau, int n, int type = 1)
	{
		double t1, t2;
		int vsize = n;
#ifdef linScheme
		vsize *= 2;
#endif

		t1 = omp_get_wtime();

		std::vector<double> w(vsize + 1), g(vsize + 2);
		int m;
		const int iterSize = 5;

		std::vector<std::vector<double>> V;
		std::vector<std::vector<double>> H;

		V.reserve(iterSize);
		H.reserve(iterSize + 1);

		V.resize(1);
		H.resize(1);

		double beta;

		std::vector <double> diag(vsize);
		for (size_t i = 0; i < n; ++i)
		{
			diag[i] = 1.0 / (2.0 * len[i]);
#ifdef linScheme
			diag[i + n] = 1.0 / (24.0 * len[i]);
#endif
		}

		V[0].resize(vsize + 1);

		std::vector<Point2D> AX(vsize);
		double timing2 = 0.0, timing3 = 0.0;
		BH.InfluenceComputation(AX, type, timing2, timing3);

		V[0] = rhs;

		SolM(V[0], len, tau, V[0], BH.pointsCopy);

		std::ofstream AXfile("AX.txt");
		AXfile.precision(16);
		for (int s = 0; s < V[0].size(); ++s)
			AXfile << V[0][s] << std::endl;
		AXfile.close();

		beta = norm(vsize + 1, V[0], '2');
		V[0] = (1 / beta) * V[0];

		double gs = beta;
		double c = 1.0;
		double s = 0.0;


		for (int j = 0; j < vsize - 1; ++j)
		{
			BH.IterativeInfluenceComputation(AX, V[j], type);

			for (int i = 0; i < n; ++i)
			{
				w[i] = (AX[i] & tau[i]) + V[j][vsize] - V[j][i] * diag[i];
#ifdef linScheme
				w[i + n] = (AX[i + n] & tau[i]) - V[j][i + n] * diag[i + n];
#endif
			}
			w[vsize] = 0.0;
			for (int i = 0; i < n; ++i)
				w[vsize] += V[j][i];

			SolM(w, len, tau, w, BH.pointsCopy);

			H.resize(j + 2);

			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			for (int i = 0; i <= j; ++i)
			{
				H[i][j] = w & V[i];
				w -= H[i][j] * V[i];
			}
			H[j + 1][j] = norm(vsize + 1, w, '2');


			if (IterRot(H, rhs, gs, c, s, j + 1, vsize))
			{
				m = j + 1;
				break;
			}

			V.push_back((1 / H[j + 1][j]) * w);
		}
		std::cout << "GMRES: " << m << std::endl;
		std::cout << "__________________________________________________" << std::endl;

		for (int i = 0; i < m + 1; i++)
			g[i] = 0.;
		g[0] = beta;

		//GivensRotations
		double oldValue;
		for (int i = 0; i < m; i++)
		{
			c = H[i][i] / sqrt(H[i][i] * H[i][i] +
				H[i + 1][i] * H[i + 1][i]);
			s = H[i + 1][i] / sqrt(H[i][i] * H[i][i] +
				H[i + 1][i] * H[i + 1][i]);

			for (int j = i + 1; j < m; j++)
			{
				oldValue = H[i][j];
				H[i][j] = c * oldValue + s * H[i + 1][j];
				H[i + 1][j] = c * H[i + 1][j] - s * oldValue;
			}


			H[i][i] = c * H[i][i] + s * H[i + 1][i];
			H[i + 1][i] = 0.;

			oldValue = g[i];
			g[i] = c * oldValue;
			g[i + 1] = -s * oldValue;
		}
		//end of GivensRotations

		std::vector<double> Y(m);
		double sum;

		// Solve HY=g
		Y[m - 1] = g[m - 1] / H[m - 1][m - 1];
		for (int k = m - 2; k >= 0; --k)
		{
			sum = 0;
			for (int s = k + 1; s < m; ++s)
				sum = sum + H[k][s] * Y[s];
			Y[k] = (g[k] - sum) / H[k][k];
		}
		// end of Solve HY=g

		for (int i = 0; i < vsize; i++)
		{
			sum = 0.0;
			for (int j = 0; j < m; j++)
				sum += V[j][i] * Y[j];

			X[i] += sum;
		}
		sum = 0.0;
		for (int j = 0; j < m; j++)
			sum += V[j][vsize] * Y[j];
		R += sum;

		t2 = omp_get_wtime();

		std::cout << (t2 - t1) - (BH.tFin - BH.tStart) << std::endl;
	}
	*/

	void SolMfile(std::vector<double>& res, const std::vector<std::vector<double>>& invM)
	{
		size_t vsize = res.size();
		std::vector<double> newres(res.size());

		for (int i = 0; i < vsize; ++i)
		{
			newres[i] = 0.0;
			for (int j = 0; j < vsize; ++j)
				newres[i] += invM[i][j] * res[j];
		}
		res = newres;
	}

	/*
	void GMRESfile(BarnesHut& BH, std::vector <double>& X, double R, std::vector<double>& rhs, std::vector<double>& len, std::vector<Point2D>& tau, int n, int type = 1)
	{
		std::vector<std::vector<double>> invM(n + 1);
		for (int i = 0; i < n + 1; ++i)
			invM[i].resize(n + 1);

		std::ifstream fM("invP-800.dat");

		for (int i = 0; i < n + 1; ++i)
			for (int j = 0; j < n + 1; ++j)
				fM >> invM[i][j];

		fM.close();

		double t1, t2;

		t1 = omp_get_wtime();

		std::vector<double> w(n + 1), g(n + 2), wold(n + 1);
		int m;
		const int iterSize = 5;

		std::vector<std::vector<double>> V;
		std::vector<std::vector<double>> H;

		V.reserve(iterSize);
		H.reserve(iterSize + 1);

		V.resize(1);
		H.resize(1);

		double beta;

		std::vector <double> diag(n);
		for (size_t i = 0; i < n / 2; ++i)
		{
			diag[i] = 1.0 / (2.0 * len[i]);
			diag[i + n / 2] = 1.0 / (24.0 * len[i]);
		}

		V[0].resize(n + 1);

		std::vector<Point2D> AX(n);

		double timing2 = 0.0, timing3 = 0.0;
		BH.InfluenceComputation(AX, type, timing2, timing3);

		for (int i = 0; i < n + 1; ++i)
		{
			V[0][i] = 0.0;
			for (int j = 0; j < n + 1; ++j)
				V[0][i] += invM[i][j] * rhs[j];
		}

		std::ofstream AXfile("AX-dir.txt");
		AXfile.precision(16);
		for (int s = 0; s < V[0].size(); ++s)
			AXfile << V[0][s] << std::endl;
		AXfile.close();

		beta = norm(n + 1, V[0], '2');
		V[0] = (1 / beta) * V[0];

		double gs = beta;
		double c = 1.0;
		double s = 0.0;


		for (int j = 0; j < n - 1; ++j)
		{
			int Vsize = (int)V.size();

			BH.IterativeInfluenceComputation(AX, V[j], type);

			for (int i = 0; i < n / 2; ++i)
			{
				wold[i] = (AX[i] & tau[i]) + V[j][n] - V[j][i] * diag[i];
				wold[i + n / 2] = (AX[i + n / 2] & tau[i]) - V[j][i + n / 2] * diag[i + n / 2];
			}

			wold[n] = 0.0;
			for (int i = 0; i < n / 2; ++i) //// 
				wold[n] += V[j][i];

			for (int i = 0; i < n + 1; ++i)
			{
				w[i] = 0.0;
				for (int j = 0; j < n + 1; ++j)
					w[i] += invM[i][j] * wold[j];
			}

			H.resize(j + 2);

			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			for (int i = 0; i <= j; ++i)
			{
				H[i][j] = w & V[i];
				w -= H[i][j] * V[i];
			}
			H[j + 1][j] = norm(n + 1, w, '2');


			if (IterRot(H, rhs, gs, c, s, j + 1, n))
			{
				m = j + 1;
				//if (j > n / 2)
				break;
			}

			V.push_back((1 / H[j + 1][j]) * w);
		}
		std::cout << "GMRES: " << m << std::endl;
		std::cout << "__________________________________________________" << std::endl;

		for (int i = 0; i < m + 1; i++)
			g[i] = 0.;
		g[0] = beta;

		//GivensRotations
		double oldValue;
		for (int i = 0; i < m; i++)
		{
			c = H[i][i] / sqrt(H[i][i] * H[i][i] +
				H[i + 1][i] * H[i + 1][i]);
			s = H[i + 1][i] / sqrt(H[i][i] * H[i][i] +
				H[i + 1][i] * H[i + 1][i]);

			for (int j = i + 1; j < m; j++)
			{
				oldValue = H[i][j];
				H[i][j] = c * oldValue + s * H[i + 1][j];
				H[i + 1][j] = c * H[i + 1][j] - s * oldValue;
			}


			H[i][i] = c * H[i][i] + s * H[i + 1][i];
			H[i + 1][i] = 0.;

			oldValue = g[i];
			g[i] = c * oldValue;
			g[i + 1] = -s * oldValue;
		}
		//end of GivensRotations

		std::vector<double> Y(m);
		double sum;

		// Solve HY=g
		Y[m - 1] = g[m - 1] / H[m - 1][m - 1];
		for (int k = m - 2; k >= 0; --k)
		{
			sum = 0;
			for (int s = k + 1; s < m; ++s)
				sum = sum + H[k][s] * Y[s];
			Y[k] = (g[k] - sum) / H[k][k];
		}
		// end of Solve HY=g

		for (int i = 0; i < n; i++)
		{
			sum = 0.0;
			for (int j = 0; j < m; j++)
				sum += V[j][i] * Y[j];

			X[i] += sum;
		}
		sum = 0.0;
		for (int j = 0; j < m; j++)
			sum += V[j][n] * Y[j];
		R += sum;

		t2 = omp_get_wtime();

		std::cout << (t2 - t1) - (BH.tFin - BH.tStart) << std::endl;
	}
	*/


	/*
	void Gauss(BarnesHut& BH, std::vector <double>& X, double R, std::vector<double>& rhs, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();

		std::vector<std::vector<double>> A;
		A.resize(2 * n + 1);
		for (int i = 0; i < 2 * n + 1; ++i)
			A[i].resize(2 * n + 1);
		double t1, t2;

		t1 = omp_get_wtime();
		BH.tree->Gauss(A, pnt);

		for (int i = 0; i < n; ++i)
		{
			A[2 * n][i] = 1.0;
			A[i][2 * n] = 1.0;
		}

		t2 = omp_get_wtime();

		std::cout << "Matrix time " << t2 - t1 << std::endl;

		std::vector<double>& b = rhs;

		std::vector <double> x(X.size() + 1);
		for (int i = 0; i < n; ++i)
			x[i] = X[i];
		x[n] = R;


		double c;

		t1 = omp_get_wtime();
		for (int k = 0; k < 2 * n; k++) // прямой ход
		{
			for (int j = k + 1; j < 2 * n + 1; j++)
			{
				c = A[j][k] / A[k][k];
				for (int i = k; i < 2 * n + 1; i++)
					A[j][i] -= c * A[k][i];
				b[j] -= c * b[k];
			}
		}

		for (int k = 2 * n; k >= 0; k--) // обратный ход
		{
			c = 0.0;
			for (int j = k + 1; j <= 2 * n; j++)
				c += A[k][j] * x[j];
			x[k] = (b[k] - c) / A[k][k];
		}

		t2 = omp_get_wtime();

		std::cout << "Gauss time " << t2 - t1 << std::endl;

		for (int i = 0; i < 2 * n; ++i)
			X[i] = x[i];
		R = x[2 * n];
	}
	*/

}//namespace BH

