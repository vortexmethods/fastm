/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
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
\version 1.3
\date 08 декабря 2022 г.
*/

#include <iostream>

namespace BH
{
	extern long long op;
	
	bool IterRot(const std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, double& c, double& s, int m, int n, double epsGMRES, int iter, bool residualShow)
	{
		bool fl;
		double hii;

		if (m == 1)
			hii = H[0][0];
		else
		{
			hii = c * H[m - 1][m - 1] - s * H[m - 2][m - 1];
			ADDOP(2);
		}
		c = hii / sqrt(hii * hii +
			H[m][m - 1] * H[m][m - 1]);
		s = H[m][m - 1] / sqrt(hii * hii +
			H[m][m - 1] * H[m][m - 1]);
		gs = -s * gs;
		ADDOP(9);

		if (residualShow)
			std::cout << "Iteration: " << iter << ", residual = " << (fabs(gs) / norm(rhs)) << std::endl;
		
		fl = ((fabs(gs) / norm(rhs)) < epsGMRES);
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
		ADDOP(10);

		for (int i = 1; i < n - 1; ++i)
		{
			a = (pnt[i].tau & pnt[i].a1);
			zn = a * alpha[i] - 1.0 / (24.0 * pnt[i].len);
			alpha[i + 1] = -(pnt[i].tau & pnt[i].c1) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			delta[i + 1] = (rhs[n + i] - a * delta[i]) / zn;
			ADDOP(12);
		}

		a = pnt[n - 2].tau & pnt[n - 2].a1;
		zn = alpha[n - 2] * a - 1.0 / (24.0 * pnt[n - 2].len);
		phi[n - 1] = -((pnt[n - 2].tau & pnt[n - 2].c1) + beta[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n + n - 2] - delta[n - 2] * a) / zn;
		ADDOP(11);

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
			ADDOP(2);
		}


		double e = (pnt[n - 1].tau & pnt[n - 1].c1);
		a = (pnt[n - 1].tau & pnt[n - 1].a1);

		yn = (rhs[n + n - 1] - e * xi[1] - a * xi[n - 1]) / (e * phi[1] + a * phi[n - 1] - 0.5 / pnt[n - 1].len);
		ADDOP(10);

		AX[n + n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
			AX[n + i] = phi[i + 1] * yn + xi[i + 1];
		ADDOP(n-1);
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
		ADDOP(11);

		for (int i = 1; i < n - 1; ++i)
		{
			a = (pnt[i].tau & pnt[i].a);
			zn = a * alpha[i] - 0.5 / pnt[i].len;
			alpha[i + 1] = -(pnt[i].tau & pnt[i].c) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			gamma[i + 1] = -(1.0 + a * gamma[i]) / zn;
			delta[i + 1] = (rhs[i] - a * delta[i]) / zn;
			ADDOP(12);
		}


		a = pnt[n - 2].tau & pnt[n - 2].a;
		zn = alpha[n - 2] * a - 0.5 / pnt[n - 2].len;
		phi[n - 1] = -((pnt[n - 2].tau & pnt[n - 2].c) + beta[n - 2] * a) / zn;
		psi[n - 1] = -(1.0 + gamma[n - 2] * a) / zn;
		xi[n - 1] = (rhs[n - 2] - delta[n - 2] * a) / zn;
		ADDOP(12);

		for (int i = n - 2; i > 0; --i)
		{
			phi[i] = alpha[i] * phi[i + 1] + beta[i];
			psi[i] = alpha[i] * psi[i + 1] + gamma[i];
			xi[i] = alpha[i] * xi[i + 1] + delta[i];
			ADDOP(3);
		}
		double e = (pnt[n - 1].tau & pnt[n - 1].c);
		a = (pnt[n - 1].tau & pnt[n - 1].a);
		lam1 = e * phi[1] + a * phi[n - 1] - 0.5 / pnt[n - 1].len;
		mu1 = e * psi[1] + a * psi[n - 1] + 1.0;
		xi1 = rhs[n - 1] - e * xi[1] - a * xi[n - 1];
		lam2 = mu2 = xi2 = 0.0;
		ADDOP(11);

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
		ADDOP(8);

#ifndef linScheme
		AX[n] = ynn;
#else
		AX[2 * n] = ynn;
#endif

		AX[n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
		{
			AX[i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];
			ADDOP(2);
		}

		//Циклическая прогонка для линейной схемы
#ifdef linScheme 
		SolCircleRun(AX, /*len, tau,*/ rhs, pnt);
#endif 
	}



	void GMRES(BarnesHut& BH, std::vector<double>& X, double R, std::vector<double>& rhs, int n, const params& prm, double& timing2, double& timing3, int& niter)
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
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (int i = 0; i < n; ++i)
		{
			diag[i] = 0.5 / BH.pointsCopyPan[i].len;
			ADDOP(1);

#ifdef linScheme
			diag[i + n] = (1.0 / 24.0) / BH.pointsCopyPan[i].len;
			ADDOP(2);
#endif
		}

#ifdef linScheme
#ifdef asympScheme
		for (int t = 0; t < prm.nAngPoints; ++t)
		{
			double cft = 0.25 * prm.mu[t] / (prm.mu[t] * prm.mu[t] - 3.0 * prm.mu[t] + 2.0);
			diag[prm.KK[t] + n] = -cft / BH.pointsCopyPan[prm.KK[t]].len;		                                       
			diag[prm.KKm[t] + n] = cft / BH.pointsCopyPan[prm.KKm[t]].len;
			ADDOP(6);
		}
#endif
#endif

		std::vector<Point2D> AX(vsize);
		double tt = 0.0;
		BH.InfluenceComputation(AX, timing2, timing3);

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

		double tPrecondStart = omp_get_wtime();
		SolM(V[0], V[0], BH.pointsCopyPan);
		double tPrecondFinish = omp_get_wtime();

		timing3 += tPrecondFinish - tPrecondStart;

		beta = norm(V[0]);
		V[0] = (1.0 / beta) * V[0];

		double gs = beta;
		double c = 1.0;
		double s = 0.0;


		for (int j = 0; j < vsize - 1; ++j)
		{
			double tt = 0.0;
			BH.IterativeInfluenceComputation(AX, V[j], timing2, timing3);

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
			
			double tPrecondStart = omp_get_wtime();
			SolM(w, w, BH.pointsCopyPan);
			double tPrecondFinish = omp_get_wtime();
			timing3 += tPrecondFinish - tPrecondStart;
			

			H.resize(j + 2);

			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			for (int i = 0; i <= j; ++i)
			{
				H[i][j] = w & V[i];
				w -= H[i][j] * V[i];
			}
			H[j + 1][j] = norm(w);

			m = j + 1;
			if (IterRot(H, rhs, gs, c, s, j + 1, vsize, BH.prm.epsGMRES, m, BH.prm.residualShow))			
				break;			

			V.push_back((1 / H[j + 1][j]) * w);
		}
		//std::cout << "GMRES: " << m <<" iterations" << std::endl;
		//printf("------------------------------------------------------------------------\n");
		niter = m;


		for (int i = 0; i < m + 1; i++)
			g[i] = 0.;
		g[0] = beta;

		//GivensRotations
		double oldValue;
		for (int i = 0; i < m; i++)
		{
			double dn = sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
			c = H[i][i] / dn;
			s = H[i + 1][i] / dn;
			ADDOP(5);

			for (int j = i + 1; j < m; j++)
			{
				oldValue = H[i][j];
				H[i][j] = c * oldValue + s * H[i + 1][j];
				H[i + 1][j] = c * H[i + 1][j] - s * oldValue;
				ADDOP(4);
			}


			H[i][i] = c * H[i][i] + s * H[i + 1][i];
			H[i + 1][i] = 0.;

			oldValue = g[i];
			g[i] = c * oldValue;
			g[i + 1] = -s * oldValue;
			ADDOP(4);
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
				sum += H[k][s] * Y[s];
			Y[k] = (g[k] - sum) / H[k][k];
			ADDOP(m-k+1);
		}
		// end of Solve HY=g

		for (int i = 0; i < vsize; i++)
		{
			sum = 0.0;
			for (int j = 0; j < m; j++)
				sum += V[j][i] * Y[j];
			ADDOP(m);
			X[i] += sum;
		}
		sum = 0.0;
		for (int j = 0; j < m; j++)
			sum += V[j][vsize] * Y[j];
		ADDOP(m);
		R += sum;

	}

}//namespace BH

