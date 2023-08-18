/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
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
\version 1.5
\date 19 июня 2024 г.
*/

#include <iostream>

namespace BH
{
	extern long long op;

	bool IterRot(std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, std::vector<double>& c, std::vector<double>& s, int m, int n, double epsGMRES, int iter, bool residualShow)
	{
		bool fl;
		double buf;

		for (int i = 0; i < m - 1; ++i)
		{
			buf = H[i][m - 1];

			H[i][m - 1] = c[i] * buf + s[i] * H[i + 1][m - 1];
			H[i + 1][m - 1] = -s[i] * buf + c[i] * H[i + 1][m - 1];
			ADDOP(4);
		}

		double zn = sqrt(H[m - 1][m - 1] * H[m - 1][m - 1] + H[m][m - 1] * H[m][m - 1]);

		c.push_back(H[m - 1][m - 1] / zn);
		s.push_back(H[m][m - 1] / zn);
		gs *= -s[m - 1];
		H[m - 1][m - 1] = c[m - 1] * H[m - 1][m - 1] + s[m - 1] * H[m][m - 1];
		H[m][m - 1] = 0.0;
		ADDOP(8);

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

		double len24 = pnt[0].len * 24.0;

		alpha[1] = (pnt[0].tau & pnt[0].c1) * len24;
		beta[1] = (pnt[0].tau & pnt[0].a1) * len24;
		delta[1] = -len24 * rhs[n + 0];
		double zn;
		ADDOP(8);

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

		yn = (rhs[n + n - 1] - e * xi[1] - a * xi[n - 1]) / (e * phi[1] + a * phi[n - 1] - 1.0 / (24.0 * pnt[n - 1].len));
		ADDOP(11);

		AX[n + n - 1] = yn;
		for (int i = n - 2; i >= 0; --i)
			AX[n + i] = phi[i + 1] * yn + xi[i + 1];
		ADDOP(n - 1);
#endif
	}


	void SolM(std::vector<double>& AX, const std::vector<double>& rhs, /*, const std::vector<PointsCopy>& pnt*/ const BarnesHut& BH, int p, int n)
	{
		//int n = (int)pnt.size();

		const std::vector<PointsCopy>& pnt = BH.pointsCopyPan[p];

		std::vector<double> alpha(n), beta(n), gamma(n), delta(n), phi(n), xi(n), psi(n);
		double yn, ynn;
		double lam1, lam2, mu1, mu2, xi1, xi2;
		double a;

		double len2 = pnt[0].len * 2.0;

		alpha[1] = (pnt[0].tau & pnt[0].c) * len2;
		beta[1] = (pnt[0].tau & pnt[0].a) * len2;
		gamma[1] = len2;
		delta[1] = -len2 * rhs[0];
		ADDOP(8);

		double zn;

		for (int i = 1; i < n - 1; ++i)
		{
			a = (pnt[i].tau & pnt[i].a);
			zn = a * alpha[i] - 0.5 / pnt[i].len;
			alpha[i + 1] = -(pnt[i].tau & pnt[i].c) / zn;
			beta[i + 1] = -a * beta[i] / zn;
			gamma[i + 1] = -(1.0 + a * gamma[i]) / zn;
			delta[i + 1] = (rhs[i] - a * delta[i]) / zn;
			ADDOP(13);
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
		for (int j = 0; j < n - 1; ++j)
		{
			lam2 += phi[j + 1];
			mu2 += psi[j + 1];
			xi2 -= xi[j + 1];
		}
		lam2 += 1.0;
		ADDOP(11);

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
			AX[i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];

		ADDOP(2 * (n - 1));

		//Циклическая прогонка для линейной схемы
#ifdef linScheme 
		SolCircleRun(AX, rhs, pnt);
#endif 
	}










#ifdef CALCSHEET
	void GMRES(BarnesHut& BH, std::vector<std::vector<double>>& X, std::vector<double>& R, std::vector<std::vector<double>>& rhs, const std::vector<int>& n, const params& prm, double& timing2, double& timing3, int& niter)
	{
		double ttGMRES1 = -omp_get_wtime();
		int totalVsize = 0;
		for (const auto& i : n)
			totalVsize += i;

#ifdef linScheme
		totalVsize *= 2;
#endif

		std::vector<double> w(totalVsize + prm.airfoilFile.size()), g(totalVsize + prm.airfoilFile.size() + 1);
		int m;
		const int iterSize = 50;

		std::vector<std::vector<double>> V;
		std::vector<std::vector<double>> H;

		V.reserve(iterSize);
		H.reserve(iterSize + 1);

		V.resize(1);
		H.resize(1);

		double beta;

		std::vector<std::vector <double>> diag(prm.airfoilFile.size());

#ifndef OLD_OMP
#pragma omp simd
#endif
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
#ifndef linScheme
			diag[p].resize(n[p]);
#else
			diag[p].resize(2 * n[p]);
#endif
			//diag[p].resize(n[p]);
			for (int i = 0; i < n[p]; ++i)
			{
				diag[p][i] = 0.5 / BH.pointsCopyPan[p][i].len;
				ADDOP(1);

#ifdef linScheme
				diag[p][i + n[p]] = (1.0 / 24.0) / BH.pointsCopyPan[p][i].len;
				ADDOP(1);
#endif
			}
			//}

#ifdef linScheme
#ifdef asympScheme
			for (int t = 0; t < prm.nAngPoints; ++t)
			{
				double cft = 0.25 * prm.mu[t] / (prm.mu[t] * prm.mu[t] - 3.0 * prm.mu[t] + 2.0);
				diag[p][prm.KK[t] + n[0]] = -cft / BH.pointsCopyPan[p][prm.KK[t]].len;
				diag[p][prm.KKm[t] + n[0]] = cft / BH.pointsCopyPan[p][prm.KKm[t]].len;
				ADDOP(6);
			}
#endif
#endif
		}

		std::vector<std::vector<std::vector<Point2D>>> AX(prm.airfoilFile.size());

		int nAX = (int)prm.airfoilFile.size();
		for (int k = 0; k < nAX; ++k)
		{
			AX[k].resize(nAX);
			for (int p = 0; p < nAX; ++p)
			{
				//AX[k][p].resize(n[k]);
#ifndef linScheme
				AX[k][p].resize(n[k]);
#else
				AX[k][p].resize(2 * n[k]);
#endif
			}
		}

		ttGMRES1 += omp_get_wtime();
		std::cout << "ttGMRES1 = " << ttGMRES1 << std::endl;


		//for (int p = 0; p < nAX * nAX; ++p)

		auto TITER = -omp_get_wtime();

		for (int i = 0; i < nAX; ++i)
			for (int j = 0; j < nAX; ++j)
			{
				BH.InfluenceComputation(AX[i][j], timing2, timing3, i, j);
			}

		TITER += omp_get_wtime();

		std::cout << "TITER_0 = " << TITER << std::endl;


		double ttGMRES2 = -omp_get_wtime();
		//Проверка, что предобуславливатель работает корректно (влияния соседних панелей рассчитаны по прямой схеме)
		for (int i = 0; i < (int)BH.pointsCopyPan.size(); ++i)
		{
			if ((BH.pointsCopyPan[0][i].a.length2() == 0.0) || (BH.pointsCopyPan[0][i].c.length2() == 0.0))
			{
				std::cout << "eps is too small!" << std::endl;
				exit(-1);
			}
		}

		ttGMRES2 += omp_get_wtime();
		std::cout << "ttGMRES2 = " << ttGMRES2 << std::endl;

		double ttGMRES3 = -omp_get_wtime();

		std::vector<std::vector<double>> residual(prm.airfoilFile.size());
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			residual[p] = rhs[p];
			for (int i = 0; i < n[p]; ++i)
			{
				residual[p][i] -= (AX[0][0][i] - prm.velInf) & BH.pointsCopyPan[p][i].tau;
				ADDOP(2);
#ifdef linScheme
				residual[p][n[p] + i] -= AX[0][0][n[p] + i] & BH.pointsCopyPan[p][i].tau;
				ADDOP(2);
#endif

				residual[p][i] += BH.pointsCopyPan[p][i].g() / (2.0 * BH.pointsCopyPan[p][i].len);
				ADDOP(2);
#ifdef linScheme
				residual[p][n[p] + i] += BH.pointsCopyPan[p][i].gamLin / (24.0 * BH.pointsCopyPan[p][i].len);
				ADDOP(2);
#endif
			}
		}


		for (int i = 0; i < prm.airfoilFile.size(); ++i)
			V[0].insert(V[0].end(), residual[i].begin(), residual[i].end());

		for (int i = 0; i < prm.airfoilFile.size(); ++i)
			V[0].push_back(0); /// суммарная гамма

		ttGMRES3 += omp_get_wtime();
		std::cout << "ttGMRES3 = " << ttGMRES3 << std::endl;



		// PRECONDITIONER
		double ttGMRES4 = -omp_get_wtime();
		//std::vector<PointsCopy> buf1;
		std::vector<double> vbuf1;
		int np = 0;

		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			//buf1.resize(n[p]);

#ifndef linScheme
			vbuf1.resize(n[p] + 1);
#else
			vbuf1.resize(2 * n[p] + 1);
#endif

#ifndef linScheme
			np += ((p == 0) ? 0 : n[p - 1]);
#else
			np += ((p == 0) ? 0 : 2 * n[p - 1]);
#endif
			for (int i = 0; i < n[p]; ++i)
			{
				//buf1[i] = BH.pointsCopyPan[p][i];
				vbuf1[i] = V[0][np + i];
#ifdef linScheme
				vbuf1[i + n[p]] = V[0][np + n[p] + i];
#endif
			}

#ifndef linScheme
			vbuf1[n[p]] = V[0][V[0].size() - (prm.airfoilFile.size() - p)];
#else
			vbuf1[2 * n[p]] = V[0][V[0].size() - (prm.airfoilFile.size() - p)];
#endif
			SolM(vbuf1, vbuf1, /*buf1*/ BH, p, n[p]);

			for (int j = np; j < np + n[p]; ++j)
			{
				V[0][j] = vbuf1[j - np];
#ifdef linScheme
				V[0][j + n[p]] = vbuf1[j - np + n[p]];
#endif
			}

			V[0][V[0].size() - (prm.airfoilFile.size() - p)] = vbuf1[vbuf1.size() - 1];
		}

		ttGMRES4 += omp_get_wtime();
		std::cout << "ttGMRES4 = " << ttGMRES4 << std::endl;

		double ttGMRES5 = -omp_get_wtime();
		beta = norm(V[0]);
		ADDOP(V[0].size() + 1);

		V[0] = (1.0 / beta) * V[0];
		ADDOP(V[0].size() + 1);

		double gs = beta;
		std::vector<double> c, s;
		c.reserve(iterSize);
		s.reserve(iterSize);
		ttGMRES5 += omp_get_wtime();
		std::cout << "ttGMRES5 = " << ttGMRES5 << std::endl;


		for (int j = 0; j < totalVsize - 1; ++j) //+ n.size()
		{
			double ttGMRES6 = -omp_get_wtime();
			BH.UpdateGams(V[j]);
			ttGMRES6 += omp_get_wtime();
			std::cout << "ttGMRES6 = " << ttGMRES6 << std::endl;

			auto TITER = -omp_get_wtime();
			for (int i = 0; i < nAX; ++i)
				for (int p = 0; p < nAX; ++p)
					BH.IterativeInfluenceComputation(AX[i][p], V[j], timing2, timing3, i, p);
			TITER += omp_get_wtime();

			std::cout << "TITER = " << TITER << std::endl;

			int cntr = 0;

			double ttGMRES7 = -omp_get_wtime();
			w.resize(0);
			w.resize(totalVsize + prm.airfoilFile.size(), 0.0);

			for (size_t pi = 0; pi < prm.airfoilFile.size(); ++pi)
			{
				for (size_t pj = 0; pj < prm.airfoilFile.size(); ++pj)
				{
#pragma omp parallel for
					for (int i = 0; i < n[pi]; ++i)
					{
						w[cntr + i] += (AX[pi][pj][i] & BH.pointsCopyPan[pi][i].tau);
						ADDOP(2);
						if (pi == pj)
						{
							w[cntr + i] += V[j][totalVsize + pi] - V[j][cntr + i] * diag[pi][i];
							ADDOP(1);
						}
#ifdef linScheme
						w[cntr + i + n[pi]] += (AX[pi][pj][i + n[pi]] & BH.pointsCopyPan[pi][i].tau);
						ADDOP(2);
						if (pi == pj)
						{
							w[cntr + i + n[pi]] += -V[j][cntr + i + n[pi]] * diag[pi][i + n[pi]];
							ADDOP(1);
						}
#endif						
					}
				}

				cntr += n[pi];
#ifdef linScheme
				cntr += n[pi];
#endif
			}

			for (int i = 0; i < n.size(); ++i)
				w[totalVsize + i] = 0.0;

			cntr = 0;
			for (int i = 0; i < n.size(); ++i)
			{
				for (int k = 0; k < n[i]; ++k)
					w[totalVsize + i] += V[j][cntr + k];
				cntr += n[i];
#ifdef linScheme
				cntr += n[i];
#endif
			}

			ttGMRES7 += omp_get_wtime();
			std::cout << "ttGMRES7 = " << ttGMRES7 << std::endl;

			// PRECONDITIONER

			double ttGMRES8 = 0.0;

			//std::vector<PointsCopy> buf2;
			std::vector<double> vbuf2;


			int np = 0;

			for (int p = 0; p < prm.airfoilFile.size(); ++p)
			{

				//buf2.resize(n[p]);

#ifndef linScheme
				vbuf2.resize(n[p] + 1);
#else
				vbuf2.resize(2 * n[p] + 1);
#endif


#ifndef linScheme
				np += ((p == 0) ? 0 : n[p - 1]);
#else
				np += ((p == 0) ? 0 : 2 * n[p - 1]);
#endif
				ttGMRES8 -= omp_get_wtime();
				for (int i = 0; i < n[p]; ++i)
				{
					//buf2[i] = BH.pointsCopyPan[p][i];
					vbuf2[i] = w[np + i];
#ifdef linScheme
					vbuf2[i + n[p]] = w[np + n[p] + i];
#endif
				}
				ttGMRES8 += omp_get_wtime();

#ifndef linScheme
				vbuf2[n[p]] = w[w.size() - (prm.airfoilFile.size() - p)];
#else
				vbuf2[2 * n[p]] = w[w.size() - (prm.airfoilFile.size() - p)];
#endif


				SolM(vbuf2, vbuf2, /*buf2*/ BH, p, n[p]);


				for (int j = np; j < np + n[p]; ++j)
				{
					w[j] = vbuf2[j - np];
#ifdef linScheme
					w[j + n[p]] = vbuf2[j - np + n[p]];
#endif
				}

				w[w.size() - (prm.airfoilFile.size() - p)] = vbuf2[vbuf2.size() - 1];
			}

			std::cout << "ttGMRES8 = " << ttGMRES8 << std::endl;



			H.resize(j + 2);


			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			double ttGMRES9 = -omp_get_wtime();
			for (int i = 0; i <= j; ++i)
			{
				double scal = 0.0;
#pragma omp parallel 
				{
#pragma omp for reduction(+:scal)
					for (int q = 0; q < (int)w.size(); ++q)
						scal += w[q] * V[i][q];

					H[i][j] = scal;// w& V[i];
#pragma omp for				
					for (int q = 0; q < (int)w.size(); ++q)
						w[q] -= scal * V[i][q];
					//w -= H[i][j] * V[i];
				}
				ADDOP(2 * V[i].size());
			}
			ttGMRES9 += omp_get_wtime();

			H[j + 1][j] = norm(w);

			ADDOP(w.size() + 1);


			//std::cout << norm(w) << std::endl;

			V.push_back((1 / H[j + 1][j]) * w);
			ADDOP(w.size() + 1);

			m = j + 1;

			if (IterRot(H, rhs[0], gs, c, s, j + 1, totalVsize, BH.prm.epsGMRES, m, BH.prm.residualShow))
				break;

			std::cout << "ttGMRES9 = " << ttGMRES9 << std::endl;

		}
		//std::cout << "GMRES: " << m <<" iterations" << std::endl;
		//printf("------------------------------------------------------------------------\n");
		niter = m;

		double ttGMRES10 = -omp_get_wtime();

		g[0] = beta;
		for (int i = 1; i < m + 1; i++)
			g[i] = 0.;


		//GivensRotations
		double oldValue;
		for (int i = 0; i < m; i++)
		{
			oldValue = g[i];
			g[i] = c[i] * oldValue;
			g[i + 1] = -s[i] * oldValue;
			ADDOP(2);
		}
		//end of GivensRotations

		std::vector<double> Y(m);
		double sum;

		// Solve HY=g
		Y[m - 1] = g[m - 1] / H[m - 1][m - 1];
		ADDOP(1);
		for (int k = m - 2; k >= 0; --k)
		{
			sum = 0;
			for (int s = k + 1; s < m; ++s)
				sum += H[k][s] * Y[s];
			Y[k] = (g[k] - sum) / H[k][k];
			ADDOP(m - k);
		}
		// end of Solve HY=g

		int cntr = 0;
		for (int p = 0; p < n.size(); p++) {
#ifndef linScheme			
			for (int i = 0; i < n[p]; i++)
#else
			for (int i = 0; i < 2 * n[p]; i++)
#endif
			{
				sum = 0.0;
				for (int j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				ADDOP(m);
				X[p][i] += sum;
			}
			cntr += n[p];
#ifdef linScheme
			cntr += n[p];
#endif
		}
		sum = 0.0;
		for (int p = 0; p < n.size(); p++) {
			for (int j = 0; j < m; j++)
				sum += V[j][totalVsize + p] * Y[j];
			ADDOP(m);
			R[p] += sum;
		}
		ttGMRES10 += omp_get_wtime();
		std::cout << "ttGMRES10 = " << ttGMRES10 << std::endl;
	}





	void GMRES_Direct(const std::vector<double>& Matr, 
		const std::vector<std::vector<double>>& Len,
		std::vector<std::vector<Point2D>>& Tau,
		std::vector<std::vector<double>>& X, std::vector<double>& R, std::vector<std::vector<double>>& rhs, const std::vector<int>& n, const params& prm, double& timing2, double& timing3, int& niter)
	{
		double ttGMRES1 = -omp_get_wtime();
		int totalVsize = 0;
		for (const auto& i : n)
			totalVsize += i;

#ifdef linScheme
		totalVsize *= 2;
#endif

		std::vector<double> w(totalVsize + prm.airfoilFile.size()), g(totalVsize + prm.airfoilFile.size() + 1);
		int m;
		const int iterSize = 50;

		std::vector<std::vector<double>> V;
		std::vector<std::vector<double>> H;

		V.reserve(iterSize);
		H.reserve(iterSize + 1);

		V.resize(1);
		H.resize(1);

		double beta;

		std::vector<std::vector <double>> diag(prm.airfoilFile.size());

#ifndef OLD_OMP
#pragma omp simd
#endif
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
#ifndef linScheme
			diag[p].resize(n[p]);
#else
			diag[p].resize(2 * n[p]);
#endif
			//diag[p].resize(n[p]);
			for (int i = 0; i < n[p]; ++i)
			{
				diag[p][i] = 0.5 / Len[p][i];
				ADDOP(1);

#ifdef linScheme
				diag[p][i + n[p]] = (1.0 / 24.0) / Len[p][i];
				ADDOP(1);
#endif
			}
			//}

#ifdef linScheme
#ifdef asympScheme
			for (int t = 0; t < prm.nAngPoints; ++t)
			{
				double cft = 0.25 * prm.mu[t] / (prm.mu[t] * prm.mu[t] - 3.0 * prm.mu[t] + 2.0);
				diag[p][prm.KK[t] + n[0]] = -cft / BH.pointsCopyPan[p][prm.KK[t]].len;
				diag[p][prm.KKm[t] + n[0]] = cft / BH.pointsCopyPan[p][prm.KKm[t]].len;
				ADDOP(6);
			}
#endif
#endif
		}

		std::vector<std::vector<std::vector<Point2D>>> AX(prm.airfoilFile.size());

		int nAX = (int)prm.airfoilFile.size();
		for (int k = 0; k < nAX; ++k)
		{
			AX[k].resize(nAX);
			for (int p = 0; p < nAX; ++p)
			{
				//AX[k][p].resize(n[k]);
#ifndef linScheme
				AX[k][p].resize(n[k]);
#else
				AX[k][p].resize(2 * n[k]);
#endif
			}
		}

		ttGMRES1 += omp_get_wtime();
		std::cout << "ttGMRES1 = " << ttGMRES1 << std::endl;


		//for (int p = 0; p < nAX * nAX; ++p)

		auto TITER = -omp_get_wtime();

		for (int i = 0; i < nAX; ++i)
			for (int j = 0; j < nAX; ++j)
			{
				//BH.InfluenceComputation(AX[i][j], timing2, timing3, i, j);
				//AX[i][j] <-- 0.0;
			}

		TITER += omp_get_wtime();

		std::cout << "TITER_0 = " << TITER << std::endl;


		double ttGMRES2 = -omp_get_wtime();
		////Проверка, что предобуславливатель работает корректно (влияния соседних панелей рассчитаны по прямой схеме)
		//for (int i = 0; i < (int)BH.pointsCopyPan.size(); ++i)
		//{
		//	if ((BH.pointsCopyPan[0][i].a.length2() == 0.0) || (BH.pointsCopyPan[0][i].c.length2() == 0.0))
		//	{
		//		std::cout << "eps is too small!" << std::endl;
		//		exit(-1);
		//	}
		//}

		ttGMRES2 += omp_get_wtime();
		std::cout << "ttGMRES2 = " << ttGMRES2 << std::endl;

		double ttGMRES3 = -omp_get_wtime();

		std::vector<std::vector<double>> residual(prm.airfoilFile.size());
		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			residual[p] = rhs[p];
			for (int i = 0; i < n[p]; ++i)
			{
				residual[p][i] -= (AX[0][0][i] - prm.velInf) & Tau[p][i];
				ADDOP(2);
#ifdef linScheme
				residual[p][n[p] + i] -= AX[0][0][n[p] + i] & Tau[p][i];
				ADDOP(2);
#endif

				residual[p][i] += /*BH.pointsCopyPan[p][i].g()*/ 0.0 / (2.0 * Len[p][i]);
				ADDOP(2);
#ifdef linScheme
				residual[p][n[p] + i] += /*BH.pointsCopyPan[p][i].gamLin*/ 0.0 / (24.0 * Len[p][i]);
				ADDOP(2);
#endif
			}
		}


		for (int i = 0; i < prm.airfoilFile.size(); ++i)
			V[0].insert(V[0].end(), residual[i].begin(), residual[i].end());

		for (int i = 0; i < prm.airfoilFile.size(); ++i)
			V[0].push_back(0); /// суммарная гамма

		ttGMRES3 += omp_get_wtime();
		std::cout << "ttGMRES3 = " << ttGMRES3 << std::endl;



		// PRECONDITIONER
		double ttGMRES4 = -omp_get_wtime();
		//std::vector<PointsCopy> buf1;
		std::vector<double> vbuf1;
		int np = 0;

		for (int p = 0; p < prm.airfoilFile.size(); ++p)
		{
			//buf1.resize(n[p]);

#ifndef linScheme
			vbuf1.resize(n[p] + 1);
#else
			vbuf1.resize(2 * n[p] + 1);
#endif

#ifndef linScheme
			np += ((p == 0) ? 0 : n[p - 1]);
#else
			np += ((p == 0) ? 0 : 2 * n[p - 1]);
#endif
			for (int i = 0; i < n[p]; ++i)
			{
				//buf1[i] = BH.pointsCopyPan[p][i];
				vbuf1[i] = V[0][np + i];
#ifdef linScheme
				vbuf1[i + n[p]] = V[0][np + n[p] + i];
#endif
			}

#ifndef linScheme
			vbuf1[n[p]] = V[0][V[0].size() - (prm.airfoilFile.size() - p)];
#else
			vbuf1[2 * n[p]] = V[0][V[0].size() - (prm.airfoilFile.size() - p)];
#endif
			//SolM(vbuf1, vbuf1, /*buf1*/ BH, p, n[p]);

			for (int j = np; j < np + n[p]; ++j)
			{
				V[0][j] = vbuf1[j - np];
#ifdef linScheme
				V[0][j + n[p]] = vbuf1[j - np + n[p]];
#endif
			}

			V[0][V[0].size() - (prm.airfoilFile.size() - p)] = vbuf1[vbuf1.size() - 1];
		}

		ttGMRES4 += omp_get_wtime();
		std::cout << "ttGMRES4 = " << ttGMRES4 << std::endl;

		double ttGMRES5 = -omp_get_wtime();
		beta = norm(V[0]);
		ADDOP(V[0].size() + 1);

		V[0] = (1.0 / beta) * V[0];
		ADDOP(V[0].size() + 1);

		double gs = beta;
		std::vector<double> c, s;
		c.reserve(iterSize);
		s.reserve(iterSize);
		ttGMRES5 += omp_get_wtime();
		std::cout << "ttGMRES5 = " << ttGMRES5 << std::endl;


		for (int j = 0; j < totalVsize - 1; ++j) //+ n.size()
		{
			double ttGMRES6 = -omp_get_wtime();
			//BH.UpdateGams(V[j]);
			ttGMRES6 += omp_get_wtime();
			std::cout << "ttGMRES6 = " << ttGMRES6 << std::endl;

			auto TITER = -omp_get_wtime();
			//for (int i = 0; i < nAX; ++i)
			//	for (int p = 0; p < nAX; ++p)
			//		BH.IterativeInfluenceComputation(AX[i][p], V[j], timing2, timing3, i, p);

			std::vector<double> mult;
			mult.resize(rhs.size(), 0.0);
			for (int ii = 0; ii < rhs.size(); ++ii)
				for (int jj = 0; jj < rhs.size(); ++jj)
					mult[ii] += Matr[ii * rhs.size() + jj] * V[j][jj];


			TITER += omp_get_wtime();

			std::cout << "TITER = " << TITER << std::endl;

			int cntr = 0;

			double ttGMRES7 = -omp_get_wtime();
			w.resize(0);
			w.resize(totalVsize + prm.airfoilFile.size(), 0.0);

			for (size_t pi = 0; pi < prm.airfoilFile.size(); ++pi)
			{
				for (size_t pj = 0; pj < prm.airfoilFile.size(); ++pj)
				{
#pragma omp parallel for
					for (int i = 0; i < n[pi]; ++i)
					{
						w[cntr + i] += (AX[pi][pj][i] & Tau[pi][i]);
						ADDOP(2);
						if (pi == pj)
						{
							w[cntr + i] += V[j][totalVsize + pi] - V[j][cntr + i] * diag[pi][i];
							ADDOP(1);
						}
#ifdef linScheme
						w[cntr + i + n[pi]] += (AX[pi][pj][i + n[pi]] & Tau[pi][i]);
						ADDOP(2);
						if (pi == pj)
						{
							w[cntr + i + n[pi]] += -V[j][cntr + i + n[pi]] * diag[pi][i + n[pi]];
							ADDOP(1);
						}
#endif						
					}
				}

				cntr += n[pi];
#ifdef linScheme
				cntr += n[pi];
#endif
			}

			for (int i = 0; i < n.size(); ++i)
				w[totalVsize + i] = 0.0;

			cntr = 0;
			for (int i = 0; i < n.size(); ++i)
			{
				for (int k = 0; k < n[i]; ++k)
					w[totalVsize + i] += V[j][cntr + k];
				cntr += n[i];
#ifdef linScheme
				cntr += n[i];
#endif
			}

			ttGMRES7 += omp_get_wtime();
			std::cout << "ttGMRES7 = " << ttGMRES7 << std::endl;

			// PRECONDITIONER

			double ttGMRES8 = 0.0;

			//std::vector<PointsCopy> buf2;
			std::vector<double> vbuf2;


			int np = 0;

			for (int p = 0; p < prm.airfoilFile.size(); ++p)
			{

				//buf2.resize(n[p]);

#ifndef linScheme
				vbuf2.resize(n[p] + 1);
#else
				vbuf2.resize(2 * n[p] + 1);
#endif


#ifndef linScheme
				np += ((p == 0) ? 0 : n[p - 1]);
#else
				np += ((p == 0) ? 0 : 2 * n[p - 1]);
#endif
				ttGMRES8 -= omp_get_wtime();
				for (int i = 0; i < n[p]; ++i)
				{
					//buf2[i] = BH.pointsCopyPan[p][i];
					vbuf2[i] = w[np + i];
#ifdef linScheme
					vbuf2[i + n[p]] = w[np + n[p] + i];
#endif
				}
				ttGMRES8 += omp_get_wtime();

#ifndef linScheme
				vbuf2[n[p]] = w[w.size() - (prm.airfoilFile.size() - p)];
#else
				vbuf2[2 * n[p]] = w[w.size() - (prm.airfoilFile.size() - p)];
#endif


				//SolM(vbuf2, vbuf2, /*buf2*/ BH, p, n[p]);


				for (int j = np; j < np + n[p]; ++j)
				{
					w[j] = vbuf2[j - np];
#ifdef linScheme
					w[j + n[p]] = vbuf2[j - np + n[p]];
#endif
				}

				w[w.size() - (prm.airfoilFile.size() - p)] = vbuf2[vbuf2.size() - 1];
			}

			std::cout << "ttGMRES8 = " << ttGMRES8 << std::endl;



			H.resize(j + 2);


			for (int i = 0; i < j + 2; ++i)
				H[i].resize(j + 1);

			double ttGMRES9 = -omp_get_wtime();
			for (int i = 0; i <= j; ++i)
			{
				double scal = 0.0;
#pragma omp parallel 
				{
#pragma omp for reduction(+:scal)
					for (int q = 0; q < (int)w.size(); ++q)
						scal += w[q] * V[i][q];

					H[i][j] = scal;// w& V[i];
#pragma omp for				
					for (int q = 0; q < (int)w.size(); ++q)
						w[q] -= scal * V[i][q];
					//w -= H[i][j] * V[i];
				}
				ADDOP(2 * V[i].size());
			}
			ttGMRES9 += omp_get_wtime();

			H[j + 1][j] = norm(w);

			ADDOP(w.size() + 1);


			//std::cout << norm(w) << std::endl;

			V.push_back((1 / H[j + 1][j]) * w);
			ADDOP(w.size() + 1);

			m = j + 1;

			if (IterRot(H, rhs[0], gs, c, s, j + 1, totalVsize, 1e-15, m, true))
				break;

			std::cout << "ttGMRES9 = " << ttGMRES9 << std::endl;

		}
		//std::cout << "GMRES: " << m <<" iterations" << std::endl;
		//printf("------------------------------------------------------------------------\n");
		niter = m;

		double ttGMRES10 = -omp_get_wtime();

		g[0] = beta;
		for (int i = 1; i < m + 1; i++)
			g[i] = 0.;


		//GivensRotations
		double oldValue;
		for (int i = 0; i < m; i++)
		{
			oldValue = g[i];
			g[i] = c[i] * oldValue;
			g[i + 1] = -s[i] * oldValue;
			ADDOP(2);
		}
		//end of GivensRotations

		std::vector<double> Y(m);
		double sum;

		// Solve HY=g
		Y[m - 1] = g[m - 1] / H[m - 1][m - 1];
		ADDOP(1);
		for (int k = m - 2; k >= 0; --k)
		{
			sum = 0;
			for (int s = k + 1; s < m; ++s)
				sum += H[k][s] * Y[s];
			Y[k] = (g[k] - sum) / H[k][k];
			ADDOP(m - k);
		}
		// end of Solve HY=g

		int cntr = 0;
		for (int p = 0; p < n.size(); p++) {
#ifndef linScheme			
			for (int i = 0; i < n[p]; i++)
#else
			for (int i = 0; i < 2 * n[p]; i++)
#endif
			{
				sum = 0.0;
				for (int j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				ADDOP(m);
				X[p][i] += sum;
			}
			cntr += n[p];
#ifdef linScheme
			cntr += n[p];
#endif
		}
		sum = 0.0;
		for (int p = 0; p < n.size(); p++) {
			for (int j = 0; j < m; j++)
				sum += V[j][totalVsize + p] * Y[j];
			ADDOP(m);
			R[p] += sum;
		}
		ttGMRES10 += omp_get_wtime();
		std::cout << "ttGMRES10 = " << ttGMRES10 << std::endl;
	}
#endif
}//namespace BH

