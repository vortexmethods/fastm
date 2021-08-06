/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: FastFourier.cpp                                                  |
| Info: Source code of FFT                                                    |
|                                                                             |
| This file is part of FFT.                                                   |
| FFT is free software: you can redistribute it and/or modify it              |
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
| along with FFT.  If not, see <http://www.gnu.org/licenses/>.                |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Реализация основного класса
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
*/


#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "FastFourier.h"


namespace FFT
{

	inline double sqr(double x) { return x * x; };

	/// Конструктор инициализации
	FastFourier::FastFourier(std::vector<Vortex2D>& wake, std::vector<Point2D>& _velo, Point2D _lD, Point2D _dim, numvector<int, 2> _nNodes, int _interpType)
		: points(wake), velo(_velo), nNodes({ _nNodes[0] + 2, _nNodes[1] + 2 }), dim(_dim), h({ _dim[0] / (_nNodes[0] - 1), _dim[1] / (_nNodes[1] - 1) }), \
		lD({ _lD[0] - _dim[0] / (_nNodes[0] - 1),_lD[1] - _dim[1] / (_nNodes[1] - 1) }), interpType(_interpType)
	{
		rU = _lD + dim + h;
		lengthOfCloseZone = 2 * sizeOfCloseZone + 1;

		Cmatr.resize(2 * (lengthOfCloseZone + 1 + 2 * interpType) * (lengthOfCloseZone + 1 + 2 * interpType) * 16);

		const int nCellsX = nNodes[0] - 1;
		const int nCellsY = nNodes[1] - 1;

		cells.resize(nCellsX * nCellsY);

		gam.resize(nNodes[0] * nNodes[1], 0.0);

#pragma omp parallel for 
		for (int i = 0; i < nCellsY; ++i)
			for (int j = 0; j < nCellsX; ++j)
			{
				int ij = i * nCellsX + j;

				cells[ij].coord = { lD[0] + j * h[0], lD[1] + i * h[1] };

				cells[ij].velCorr.toZero(Point2D({ 0.0, 0.0 }));

				cells[ij].velDerX.toZero(Point2D({ 0.0, 0.0 }));
				cells[ij].velDerY.toZero(Point2D({ 0.0, 0.0 }));

				cells[ij].closeZone.reserve(lengthOfCloseZone * lengthOfCloseZone);
				cells[ij].numOfCloseZone.reserve(lengthOfCloseZone * lengthOfCloseZone);

				cells[ij].index.reserve((int)(points.size() / (nCellsX * nCellsY)));

				int num = 0;


				for (int k = 0; k < lengthOfCloseZone; ++k)
				{
					int bou1 = (i + k - sizeOfCloseZone);
					for (int m = 0; m < lengthOfCloseZone; m++)
					{
						int bou2 = (j + m - sizeOfCloseZone);

						int ind = bou1 * nCellsX + bou2;

						if ((bou1 >= 0) && (bou1 < nCellsY) && (bou2 >= 0) && (bou2 < nCellsX))
						{
							cells[ij].closeZone.push_back(&(cells[ind]));
							cells[ij].numOfCloseZone.push_back(num);
						}
						num++;
					}// for m
				}// for k
			}// for j
	}//конструктор

	/// Одномерный оператор сглаживания
	double FastFourier::M4(double x)
	{
		double x2 = x * x;
		double abs_x = fabs(x);
		if (abs_x <= 1.0)
			return 1.0 - 2.5 * x2 + 1.5 * x2 * abs_x;

		else if (abs_x <= 2.0)
			return 0.5 * (2.0 - abs_x) * (2.0 - abs_x) * (1.0 - abs_x);

		else
			return 0.0;
	}//M4(...)

	/// Двумерный оператор сглаживания 
	double FastFourier::W(const Point2D& r)
	{
		return M4(r[0] / h[0]) * M4(r[1] / h[1]);
	}//W(...)

	void FastFourier::InitializationCells()
	{
		int i, j;
		int cnt = 0;

		for (auto& it : points)
		{
			i = (int)((it.r()[0] - lD[0]) / h[0]);
			j = (int)((it.r()[1] - lD[1]) / h[1]);
			cells[j * (nNodes[0] - 1) + i].index.push_back(cnt);
			cnt++;
		}
	}

	/// Вычисление матрицы коррекции С
	void FastFourier::ReadMatrC()
	{
		std::stringstream ss;
		ss << "../../test/CMatr/CMatr";
		ss << sizeOfCloseZone + interpType;
		ss << "_";
		ss << nNodes[1];
		ss << ".dat";

		std::ifstream matrfile(ss.str());
		if (!matrfile.good())
		{
			printf("!!! Error reading file with mesh !!!");
			exit(99);
		}

		int Csize = 2 * (lengthOfCloseZone + 1 + 2 * interpType) * (lengthOfCloseZone + 1 + 2 * interpType);

		for (int i = 0; i < Csize; i++)
			for (int j = 0; j < 16; j++)
				matrfile >> Cmatr[i * 16 + j];

		matrfile.close();
	}//ReadMatrC()

	/// Вычисление циркуляций в 16 узлах сетки ("размазывание" влияния ячейки q)
	numvector<double, 16> FastFourier::CalcGam16(const cell& q)
	{
		numvector<double, 16> gam16;
		gam16.toZero();

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				Point2D qq = q.coord + Point2D({ h[0] * (j - 1), h[1] * (i - 1) });

				for (size_t qi : q.index)
				{
					const Vortex2D& itDop = points[qi];
					Point2D vec = itDop.r() - qq;
					gam16[i * 4 + j] += itDop.g() * W(vec);
				}
			}
		return gam16;
	}//CalcGam16(...)

	/// Вычисление вектора циркуляций в узлах сетки с помощью оператора сглаживания
	void FastFourier::FillGam()
	{
		const int nCellsX = nNodes[0] - 1;
		const int nCellsY = nNodes[1] - 1;

		gam.resize(nNodes[0] * nNodes[1], 0.0);

#pragma omp parallel for
		for (int i = 1; i < nCellsY - 1; i++)
			for (int j = 1; j < nCellsX - 1; j++)
			{
				cell& q = cells[i * nCellsX + j];
				q.gam16 = CalcGam16(q);

				for (int k = 0; k < 4; k++)
					for (int m = 0; m < 4; m++)
#pragma omp atomic
						gam[(i - 1 + k) * nNodes[0] + j - 1 + m] += q.gam16[4 * k + m];

			}//for i,j
	}//FillGam(...)


	/// Вычисление коррекционных скоростей, заполнение velCorr для каждой ячейки
	void FastFourier::CalcCorrectionVel()
	{
		std::vector<Point2D> dopVel;
		dopVel.resize((lengthOfCloseZone + 1 + 2 * interpType) * (lengthOfCloseZone + 1 + 2 * interpType), { 0.0, 0.0 });

		numvector<double, 16> gam16;

		int index;
		const int nCellsX = nNodes[0] - 1;
		const int nCellsY = nNodes[1] - 1;
		const int nC = nCellsX * nCellsY;

		std::vector<Point2D> dopVelCorr, dopVelDerX, dopVelDerY;

		size_t thnum = 0;

#pragma omp parallel shared(dopVelCorr, dopVelDerX, dopVelDerY, thnum) private(gam16, index) firstprivate(dopVel) //num_threads(1)
		{
			size_t ompth = omp_get_thread_num();
			thnum = omp_get_num_threads();

#pragma omp single
			dopVelCorr.resize(thnum * nC * 4, { 0.0, 0.0 });

			if (interpType == 1)
			{
#pragma omp single
				dopVelDerX.resize(thnum * nC * 4, { 0.0, 0.0 });

#pragma omp single
				dopVelDerY.resize(thnum * nC * 4, { 0.0, 0.0 });
			}


#pragma omp for 
			for (int i = 1; i < nCellsY - 1; ++i)
				for (int j = 1; j < nCellsX - 1; ++j)
				{
					index = i * nCellsX + j;
					const cell& q = cells[index];
					gam16 = q.gam16;

					for (size_t s = 0; s < dopVel.size(); ++s)
					{
						double tx = 0, ty = 0;

						for (int m = 0; m < 16; ++m)
							tx += Cmatr[32 * s + m] * gam16[m];

						for (int m = 0; m < 16; ++m)
							ty += Cmatr[32 * s + 16 + m] * gam16[m];

						dopVel[s][0] = tx;
						dopVel[s][1] = ty;
					}

					int ind, ind2, ind3, ind4, ii;

					int sss = 0; // индекс по вектору номеров ближней зоны
					for (int k = 0; k < lengthOfCloseZone; ++k)
						for (int m = 0; m < lengthOfCloseZone; ++m)
						{
							ind = (k + interpType) * (lengthOfCloseZone + 1 + 2 * interpType) + (m + interpType);
							ind2 = ind + (lengthOfCloseZone + 1 + 2 * interpType); //строка выше
							ind3 = ind - (lengthOfCloseZone + 1 + 2 * interpType); //строка ниже 
							ind4 = ind2 + (lengthOfCloseZone + 1 + 2 * interpType); //2 строки выше

							if (sss < int(q.numOfCloseZone.size()))
							{
								int indCZ = k * lengthOfCloseZone + m;
								if (indCZ == q.numOfCloseZone[sss])
								{
									int locind = (i - sizeOfCloseZone + k) * nCellsX + (j - sizeOfCloseZone + m);

									ii = (int)((ompth * nC + locind) * 4);

									dopVelCorr[ii + 0] += dopVel[ind];
									dopVelCorr[ii + 1] += dopVel[ind + 1];
									dopVelCorr[ii + 2] += dopVel[ind2];
									dopVelCorr[ii + 3] += dopVel[ind2 + 1];

									if (interpType == 1)
									{
										dopVelDerX[ii + 0] += dopVel[ind + 1] - dopVel[ind - 1];
										dopVelDerX[ii + 1] += dopVel[ind + 2] - dopVel[ind];
										dopVelDerX[ii + 2] += dopVel[ind2 + 1] - dopVel[ind2 - 1];
										dopVelDerX[ii + 3] += dopVel[ind2 + 2] - dopVel[ind2];

										dopVelDerY[ii + 0] += dopVel[ind2] - dopVel[ind3];
										dopVelDerY[ii + 1] += dopVel[ind2 + 1] - dopVel[ind3 + 1];
										dopVelDerY[ii + 2] += dopVel[ind4] - dopVel[ind];
										dopVelDerY[ii + 3] += dopVel[ind4 + 1] - dopVel[ind + 1];
									}
									++sss;
								}
							}//if sss
						}//for k,m
				}//for i,j
		}//parallel section

#pragma omp parallel for 
		for (int j = 0; j < 4; j++)
		{
			int ind;
			for (size_t i = 0; i < cells.size(); i++)
				for (size_t k = 0; k < thnum; k++)
				{
					ind = (int)((k * nC + i) * 4 + j);
					cells[i].velCorr[j] -= dopVelCorr[ind];
					if (interpType == 1)
					{
						cells[i].velDerX[j] -= 0.5 / h[0] * dopVelDerX[ind];
						cells[i].velDerY[j] -= 0.5 / h[1] * dopVelDerY[ind];
					}
				}
		}
	}//CalcCorrectionVel(...)


	/// Расчет скоростей ВЭ внутри q ячейки по закону Био --- Савара от ближней зоны
	inline void FastFourier::CalcVeloByBiotSavartForCell(const cell& q)
	{
		//Локальные переменные для цикла
		Point2D velI;
		Point2D tempVel;
		double dst2eps;
		const double cft = 0.5 / PI;

		//сама q ячейка входит в вектор CloseZone
		for (size_t k = 0; k < q.closeZone.size(); ++k)
		{
			//ячейка, от которой считаем влияние 
			const cell& p = *(q.closeZone[k]);

			for (int i : q.index)
			{
				velI.toZero();

				const Point2D& posI = points[i].r();

				for (int j : p.index)
				{
					const Point2D& posJ = points[j].r();
					const double& gamJ = points[j].g();

					dst2eps = std::max(dist2(posI, posJ), eps2);

					tempVel = { (-posI[1] + posJ[1]), (posI[0] - posJ[0]) };
					velI += tempVel * (gamJ / dst2eps);
				}//for j

				velI *= cft;

				velo[i] += velI;
			}//for i 
		}//for k

	}//CalcVeloByBiotSavartForCell(...)


	/// Функция ядра Био --- Савара
	Point2D FastFourier::Skos(const Point2D& R, const Point2D& X) const
	{
		double dst2 = std::max(dist2(R, X), eps2);
		Point2D res = { -(R[1] - X[1]), (R[0] - X[0]) };
		res *= IDPI / dst2;
		return res;
	}//Skos(...)


	/// решение уравнение Пуассона, заполнение вектора velPuasson
	void FastFourier::PuassonSolver(double& tExclude)
	{
		double t1, t2;
		tExclude = 0.0;
		Eigen::MatrixXd omega, bkx, bky;
		omega.setZero(2 * nNodes[0], 2 * nNodes[1]);
		t1 = omp_get_wtime();
		bkx.setZero(2 * nNodes[0], 2 * nNodes[1]);
		bky.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int i = 0; i < nNodes[0]; ++i)
			for (int j = 0; j < nNodes[1]; ++j)
			{
				omega(j, i) = gam[i * nNodes[1] + j];
				Point2D skos = Skos(Point2D({ i * h[0], j * h[1] }), { 0.0, 0.0 });
				bkx(i, j) = skos[0];
				bky(i, j) = skos[1];
			}

		for (int i = 0; i < nNodes[0] - 1; i++)
			for (int j = 0; j < nNodes[1]; j++)
			{
				bkx(2 * nNodes[0] - 1 - i, j) = bkx(i + 1, j);
				bky(2 * nNodes[0] - 1 - i, j) = -bky(i + 1, j);
			}

		for (int i = 0; i < nNodes[0]; i++)
			for (int j = 0; j < nNodes[1] - 1; j++)
			{
				bkx(i, 2 * nNodes[1] - 1 - j) = -bkx(i, j + 1);
				bky(i, 2 * nNodes[1] - 1 - j) = bky(i, j + 1);
			}

		for (int i = 0; i < nNodes[0] - 1; i++)
			for (int j = 0; j < nNodes[1] - 1; j++)
			{
				bkx(2 * nNodes[0] - 1 - i, 2 * nNodes[1] - 1 - j) = -bkx(i + 1, j + 1);
				bky(2 * nNodes[0] - 1 - i, 2 * nNodes[1] - 1 - j) = -bky(i + 1, j + 1);
			}
		t2 = omp_get_wtime();
		tExclude += t2 - t1;

		velPuasson.resize(nNodes[0] * nNodes[1], { 0.0, 0.0 });

		//FFT

		Eigen::FFT<double> fft;

		Eigen::MatrixXcd fOmega;
		fOmega.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int k = 0; k < omega.rows(); k++) {
			Eigen::VectorXcd tmpOmega(2 * nNodes[0]);
			fft.fwd(tmpOmega, omega.row(k));
			fOmega.row(k) = tmpOmega;
		}

		for (int k = 0; k < omega.cols(); k++) {
			Eigen::VectorXcd tmpOmega(2 * nNodes[1]);
			fft.fwd(tmpOmega, fOmega.col(k));
			fOmega.col(k) = tmpOmega;
		}

		t1 = omp_get_wtime();
		Eigen::MatrixXcd fBkx;
		fBkx.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int k = 0; k < bkx.rows(); k++) {
			Eigen::VectorXcd tmpBkx(2 * nNodes[0]);
			fft.fwd(tmpBkx, bkx.row(k));
			fBkx.row(k) = tmpBkx;
		}

		for (int k = 0; k < bkx.cols(); k++) {
			Eigen::VectorXcd tmpBkx(2 * nNodes[1]);
			fft.fwd(tmpBkx, fBkx.col(k));
			fBkx.col(k) = tmpBkx;
		}

		Eigen::MatrixXcd fBky;
		fBky.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int k = 0; k < bky.rows(); k++) {
			Eigen::VectorXcd tmpBky(2 * nNodes[0]);
			fft.fwd(tmpBky, bky.row(k));
			fBky.row(k) = tmpBky;
		}

		for (int k = 0; k < bky.cols(); k++) {
			Eigen::VectorXcd tmpBky(2 * nNodes[1]);
			fft.fwd(tmpBky, fBky.col(k));
			fBky.col(k) = tmpBky;
		}
		t2 = omp_get_wtime();
		tExclude += t2 - t1;

		Eigen::MatrixXcd fUx, fUy;
		fUx.setZero(2 * nNodes[0], 2 * nNodes[1]);
		fUy.setZero(2 * nNodes[0], 2 * nNodes[1]);

#pragma omp parallel for 
		for (int i = 0; i < 2 * nNodes[0]; ++i)
			for (int j = 0; j < 2 * nNodes[1]; ++j)
			{
				fUx(i, j) = fBkx(i, j) * fOmega(i, j);
				fUy(i, j) = fBky(i, j) * fOmega(i, j);
			}

		Eigen::MatrixXcd ux;
		ux.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int k = 0; k < fUx.rows(); k++) {
			Eigen::VectorXcd tmpUx(2 * nNodes[0]);
			fft.inv(tmpUx, fUx.row(k));
			ux.row(k) = tmpUx;
		}

		for (int k = 0; k < fUx.cols(); k++) {
			Eigen::VectorXcd tmpUx(2 * nNodes[1]);
			fft.inv(tmpUx, ux.col(k));
			ux.col(k) = tmpUx;
		}

		Eigen::MatrixXcd uy;
		uy.setZero(2 * nNodes[0], 2 * nNodes[1]);

		for (int k = 0; k < fUy.rows(); k++) {
			Eigen::VectorXcd tmpUy(2 * nNodes[0]);
			fft.inv(tmpUy, fUy.row(k));
			uy.row(k) = tmpUy;
		}

		for (int k = 0; k < fUy.cols(); k++) {
			Eigen::VectorXcd tmpUy(2 * nNodes[1]);
			fft.inv(tmpUy, uy.col(k));
			uy.col(k) = tmpUy;
		}

		for (int i = 0; i < nNodes[1]; ++i)
			for (int j = 0; j < nNodes[0]; ++j)
				velPuasson[i * nNodes[0] + j] = { ux(i, j).real(), uy(i, j).real() };
	}//PuassonSolver()


	/// Линейная интерполяция на ВЭ в q-й ячейке
	void FastFourier::InterpolationFromCell(cell& q)
	{
		double idh = 1.0 / (h[0] * h[1]);

		for (int i : q.index)
		{
			Point2D locXY = points[i].r() - q.coord;
			velo[i] += idh * (q.velCorr[0] * (h[0] - locXY[0]) * (h[1] - locXY[1]) + \
				q.velCorr[2] * (h[0] - locXY[0]) * locXY[1] + \
				q.velCorr[1] * locXY[0] * (h[1] - locXY[1]) + \
				q.velCorr[3] * locXY[0] * locXY[1]);
		}

	}//InterpolationFromCell(...)

	/// Интерполяция Эрмита на ВЭ в q-й ячейке
	void FastFourier::HermiteInterpolationFromCell(cell& q)
	{
		double idh = 1.0 / (h[0] * h[1]);

		// Отклонения от левого нижнего узла ячейки
		double dx, dy;

		double hx[4];
		double hy[4];
		double H[12];

		for (int i : q.index)
		{
			dx = (points[i].r()[0] - q.coord[0]) / h[0];
			dy = (points[i].r()[1] - q.coord[1]) / h[1];

			hx[0] = sqr(1.0 - dx) * (1.0 + 2.0 * dx);
			hx[1] = sqr(dx) * (3.0 - 2.0 * dx);
			hx[2] = dx * sqr(1.0 - dx);
			hx[3] = sqr(dx) * (dx - 1.0);

			hy[0] = sqr(1.0 - dy) * (1.0 + 2.0 * dy);
			hy[1] = sqr(dy) * (3.0 - 2.0 * dy);
			hy[2] = dy * sqr(1.0 - dy);
			hy[3] = sqr(dy) * (dy - 1.0);

			H[0] = hx[0] * hy[0];
			H[1] = hx[1] * hy[0];
			H[2] = hx[0] * hy[1];
			H[3] = hx[1] * hy[1];

			H[4] = hx[2] * hy[0];
			H[5] = hx[3] * hy[0];
			H[6] = hx[2] * hy[1];
			H[7] = hx[3] * hy[1];

			H[8] = hx[0] * hy[2];
			H[9] = hx[1] * hy[2];
			H[10] = hx[0] * hy[3];
			H[11] = hx[1] * hy[3];

			velo[i] += q.velCorr[0] * H[0] + q.velCorr[1] * H[1] + q.velCorr[2] * H[2] + q.velCorr[3] * H[3] + \
				h[0] * (q.velDerX[0] * H[4] + q.velDerX[1] * H[5] + q.velDerX[2] * H[6] + q.velDerX[3] * H[7]) + \
				h[1] * (q.velDerY[0] * H[8] + q.velDerY[1] * H[9] + q.velDerY[2] * H[10] + q.velDerY[3] * H[11]);
		}

	}//InterpolationFromCell(...)

	/// Расчет итоговых скоростей
	void FastFourier::CalcTotalVelo()
	{
		numvector<Point2D, 4> velP, velDerX, velDerY;

		//в пустых не считаем
#pragma omp parallel for private(velP, velDerX, velDerY) 
		for (int i = 1; i < nNodes[1] - 2; i++)
			for (int j = 1; j < nNodes[0] - 2; j++)
			{
				int indC = i * (nNodes[0] - 1) + j;

				velP[0] = velPuasson[j * nNodes[1] + i];
				velP[2] = velPuasson[j * nNodes[1] + i + 1];
				velP[1] = velPuasson[(j + 1) * nNodes[1] + i];
				velP[3] = velPuasson[(j + 1) * nNodes[1] + i + 1];

				cells[indC].velCorr += velP;

				if (interpType == 1)
				{
					velDerX[0] = velPuasson[(j + 1) * nNodes[1] + i] - velPuasson[(j - 1) * nNodes[1] + i];
					velDerX[1] = velPuasson[(j + 2) * nNodes[1] + i] - velPuasson[j * nNodes[1] + i];
					velDerX[2] = velPuasson[(j + 1) * nNodes[1] + i + 1] - velPuasson[(j - 1) * nNodes[1] + i + 1];
					velDerX[3] = velPuasson[(j + 2) * nNodes[1] + i + 1] - velPuasson[j * nNodes[1] + i + 1];

					velDerY[0] = velPuasson[j * nNodes[1] + i + 1] - velPuasson[j * nNodes[1] + i - 1];
					velDerY[1] = velPuasson[(j + 1) * nNodes[1] + i + 1] - velPuasson[(j + 1) * nNodes[1] + i - 1];
					velDerY[2] = velPuasson[j * nNodes[1] + i + 2] - velPuasson[j * nNodes[1] + i];
					velDerY[3] = velPuasson[(j + 1) * nNodes[1] + i + 2] - velPuasson[(j + 1) * nNodes[1] + i];

					cells[indC].velDerX += 0.5 / h[0] * velDerX;
					cells[indC].velDerY += 0.5 / h[1] * velDerY;

					HermiteInterpolationFromCell(cells[indC]);
				}
				else 
					InterpolationFromCell(cells[indC]);

				CalcVeloByBiotSavartForCell(cells[indC]);
			}
	}//CalcTotalVelo()

}//namespace FFT