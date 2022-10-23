/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.2    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/10/22     |
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
\version 1.2
\date 22 октября 2022 г.
*/

#include <iostream>
#include "BarnesHut.h"

namespace BH
{
	
	
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

		if ((fabs(gs) / norm(n, rhs, '2')) < epsGMRES) fl = true;
		else fl = false;
		return fl;
	}

	void SolCircleRun(std::vector<double>& AX, const std::vector<double>& rhs, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();
#ifdef linScheme
		std::vector<double> alpha(n), beta(n), delta(n), phi(n), xi(n);
		double yn, a;

		alpha[1] = (pnt[0].tau & pnt[0].c1) * pnt[0].len * 24.0;
		beta[1] = -(pnt[n - 1].tau & pnt[0].a1) * pnt[0].len * 24.0; ///почему тут минус?
		delta[1] = -24.0 * pnt[0].len * rhs[n + 0];
		double zn;

		for (int i = 1; i < n - 1; ++i)
		{
			a = (pnt[i].tau & pnt[i].a1);
			zn = a * alpha[i] - 1.0 / (24.0 * pnt[i].len);
			alpha[i + 1] = -(pnt[i].tau & pnt[i].c1) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			delta[i + 1] = (rhs[n + i] - a * delta[i]) / zn;
		}

		a = pnt[n - 2].tau & pnt[n - 2].a1;
		zn = alpha[n - 2] * a - 1.0 / (24.0 * pnt[n - 2].len);
		phi[n - 1] = -((pnt[n - 2].tau & pnt[n - 2].c1) + beta[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n + n - 2] - delta[n - 2] * a) / zn;

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
		}

		double e = (pnt[n - 1].tau & pnt[n - 1].c1);
		a = (pnt[n - 1].tau & pnt[n - 1].a1);

		yn = (rhs[n + n - 1] - e * xi[1] - a * xi[n - 1]) / (e * phi[1] + a * phi[n - 1] - 0.5 / pnt[n - 1].len);

		AX[n + n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
			AX[n + i] = phi[i + 1] * yn + xi[i + 1];
#endif
	}


	void SolM(std::vector<double>& AX, const std::vector<double>& rhs, const std::vector<PointsCopy>& pnt)
	{
		int n = (int)pnt.size();

		std::vector<double> alpha(n), beta(n), gamma(n), delta(n), phi(n), xi(n), psi(n);
		double yn, ynn;
		double lam1, lam2, mu1, mu2, xi1, xi2;
		double a;

		alpha[1] = (pnt[0].tau & pnt[0].c) * pnt[0].len * 2.0;
		beta[1] = -(pnt[n - 1].tau & pnt[0].a) * pnt[0].len * 2.0; ///почему тут минус?
		gamma[1] = 2.0 * pnt[0].len;
		delta[1] = -2.0 * pnt[0].len * rhs[0];
		double zn;

		for (int i = 1; i < n - 1; ++i)
		{
			a = (pnt[i].tau & pnt[i].a);
			zn = a * alpha[i] - 0.5 / pnt[i].len;
			alpha[i + 1] = -(pnt[i].tau & pnt[i].c) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			gamma[i + 1] = -(1.0 + a * gamma[i]) / zn;
			delta[i + 1] = (rhs[i] - a * delta[i]) / zn;
		}

		a = pnt[n - 2].tau & pnt[n - 2].a;
		zn = alpha[n - 2] * a - 0.5 / pnt[n - 2].len;
		phi[n - 1] = -((pnt[n - 2].tau & pnt[n - 2].c) + beta[n - 2] * a) / zn;
		psi[n - 1] = -(1.0 + gamma[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n - 2] - delta[n - 2] * a) / zn;

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			psi[i] = alpha[i] * psi[i + 1] + gamma[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
		}
		double e = (pnt[n - 1].tau & pnt[n - 1].c);
		a = (pnt[n - 1].tau & pnt[n - 1].a);
		lam1 = e * phi[1] + a * phi[n - 1] - 0.5 / pnt[n - 1].len;
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
		SolCircleRun(AX, /*len, tau,*/ rhs, pnt);
#endif 
	}


	
	void GMRES(BarnesHut& BH, std::vector<double>& X, double R, std::vector<double>& rhs, int n)
	{
		int vsize = n;
#ifdef linScheme
		vsize *= 2;
#endif

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
		for (int i = 0; i < n; ++i)
		{
			diag[i] = 0.5 / BH.pointsCopyPan[i].len;
#ifdef linScheme
			diag[i + n] = (1.0 / 24.0) / BH.pointsCopyPan[i].len;
#endif
		}

#ifdef linScheme
#ifdef asympScheme
		for (int t = 0; t < nAngPoints; ++t)
		{
			double cft = 0.25 * mu[t] / (mu[t] * mu[t] - 3.0 * mu[t] + 2.0);
			diag[KK[t] + n] = -cft / BH.pointsCopyPan[KK[t]].len;
		                                          //len[KK[t]];
			diag[KKm[t] + n] = cft / BH.pointsCopyPan[KKm[t]].len;
												//len[KKm[t]];
		}
#endif
#endif

		std::vector<Point2D> AX(vsize);
		double timing2 = 0.0, timing3 = 0.0;
		BH.InfluenceComputation(AX, timing2, timing3);
		//std::cout << "time2 = " << timing2 << ", time3 = " << timing3 << std::endl;


		//Проверка, что предобуславливатель работает корректно (влияния соседних панелей рассчитаны по прямой схеме)
		for (int i = 0; i < (int)BH.pointsCopyPan.size(); ++i)
		{
			if ((BH.pointsCopyPan[i].a.length2() == 0.0) || (BH.pointsCopyPan[i].c.length2() == 0.0))
			{
				std::cout << "eps is too small!" << std::endl;
				exit(-1);
			}
		}


		V[0] = rhs;
		V[0].push_back(0); /// суммарная гамма

		SolM(V[0], V[0], BH.pointsCopyPan);

		beta = norm(vsize + 1, V[0], '2');
		V[0] = (1.0 / beta) * V[0];

		double gs = beta;
		double c = 1.0;
		double s = 0.0;


		for (int j = 0; j < vsize - 1; ++j)
		{
			//std::cout << "j = " << j << std::endl;
			timing2 = 0.0; timing3=0.0;
			BH.IterativeInfluenceComputation(AX, V[j], timing2, timing3);

			//std::cout << "time2 = " << timing2 << ", time3 = " << timing3 << std::endl;

			for (int i = 0; i < n; ++i)
			{
				w[i] = (AX[i] & BH.pointsCopyPan[i].tau) + V[j][vsize] - V[j][i] * diag[i];
#ifdef linScheme
				w[i + n] = (AX[i + n] & BH.pointsCopyPan[i].tau) - V[j][i + n] * diag[i + n];
#endif
			}
			w[vsize] = 0.0;
			for (int i = 0; i < n; ++i)
				w[vsize] += V[j][i];

			SolM(w, w, BH.pointsCopyPan);

			H.resize(j + 2);

			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			for (int i = 0; i <= j; ++i)
			{
				H[i][j] = w & V[i];
				w -= H[i][j] * V[i];
			}
			H[j + 1][j] = norm(vsize + 1, w, '2');

			m = j + 1;
			if (IterRot(H, rhs, gs, c, s, j + 1, vsize))			
				break;			

			V.push_back((1 / H[j + 1][j]) * w);
		}
		std::cout << "GMRES: " << m <<" iterations" << std::endl;
		printf("------------------------------------------------------------------------\n");

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

	}

}//namespace BH

