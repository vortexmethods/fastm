/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: FastFourier.h                                                    |
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
\brief �������� ��������� ������
\author ���������� ���� ��������������
\author ������ ������� ��������
\version 1.0
\date 05 ������� 2021 �.
*/

#ifndef FASTFOURIER_H_
#define FASTFOURIER_H_

#include <vector>

//#define EIGEN_FFTW_DEFAULT
#include "Eigen/Eigen"
#include "unsupported/Eigen/FFT"

#include "Params.h"
#include "Vortex2D.h"


namespace FFT
{

	const double PI = 3.1415926535897932;
	const double IDPI = 0.15915494309189534;

	/// �������� ��������� ����� ������ �����
	struct cell
	{
		/// ���������� ������ ������� ���� ������
		Point2D coord;

		/// ������ �� ���������� �� ������ ������� ����
		std::vector<cell*> closeZone;

		/// ������ ������� �������� ����� ������� ���� � ����� �������
		std::vector<int> numOfCloseZone;

		/// ������ ������������� ��������� � ����� ������ 
		numvector<Point2D, 4> velCorr;

		/// ������ ����������� �������� � ����� ������ 
		numvector<Point2D, 4> velDerX, velDerY;

		/// ������� ������ ������ ������ � ����� �������
		std::vector<int> index;

		/// ������� ������, ���������� � "�������" ������� ���� (������ � �����)
		std::vector<int> indexRight, indexLeft;

		/// ������ ���������� � 9-�� �������� �������, �� ������� "�������������" ������� ������
		numvector<double, 16> gam16;

	};



	class FastFourier
	{

	private:
		/// ����� "������� ����" (� �������)
		int lengthOfCloseZone;

		/// ������ ��������� �������� ������, ������������ ��� ���������� ������� ���������
		std::vector<Point2D> posForC;

		/// ������� ��������� �
		//std::vector<numvector<double, 16>> Cmatr; 
		std::vector<double> Cmatr;

		/// ����� ����� ����� (Nx - �� �����������, Ny - �� ��������� )
		const numvector<int, 2> nNodes;

		/// ���� �����
		const Point2D h;

		/// ���������� ������ ������� ���� ������� (const)
		const Point2D lD;

		/// ���������� ������� �������� ���� ������� (const)
		Point2D rU;

		/// ������� �������
		// dim[0] --- �� ����������� 
		// dim[1] --- �� ���������
		const Point2D dim;

		/// ��� ������������ � ����� ����� �� �� (0 --- ����������, 1 --- ������������ ������)
		const int interpType;

	public:
		/// ������ �� ������
		std::vector<Vortex2D>& points;

		/// ������ �� ���������
		std::vector<Point2D>& velo;

		/// ������ �� ����� �����
		std::vector<cell> cells;

		/// ������ ���������� � ����� ����� (���� ��� ������� �������� MPI)
		std::vector<double> gam;

		/// ������ ���������, ���������� �������� ��������� �������� �� ���� ����� 
		std::vector<Point2D> velPuasson;

		std::vector<double> bkx;
		std::vector<double> bky;
		std::vector<double> omega;
		std::vector<double> hUx;
		std::vector<double> hUy;


		/// ����������� �������������
		FastFourier(std::vector<Vortex2D>& wake, std::vector<Point2D>& _velo, Point2D _lD, Point2D _dim, numvector<int, 2> _nNodes, int _interpType);

		/// ���������� �������� ����������� 
		double M4(double x);

		/// ��������� �������� �����������
		double W(const Point2D& r);

		/// ���������� ������� ���������� � ����� ����� � ������� ��������� �����������
		void FillGam();

		/// ������������� ����� �����
		void InitializationCells();

		/// ���������� ������� ��������� � 
		// ���������� �� �����
		void ReadMatrC();

		/// ���������� ������������� ��������� � �����������, ���������� Velcorr, VelDerX � VelDerY ��� ������ ������
		void CalcCorrectionVel();

		/// ������ ��������� �� ������ q ������ �� ������ ��� --- ������ �� ������� ����
		inline void CalcVeloByBiotSavartForCell(const cell& q);

		/// ������� ���� ��� --- ������
		Point2D Skos(const Point2D& R, const Point2D& X) const;

		/// ������� ��������� ��������, ���������� ������� VelPuasson
		void PuassonSolver(double& tExclude);

		/// �������� ������������ � ����� �� �� � q-� ������
		void InterpolationFromCell(cell& q);

		/// ������������ ������ � ����� �� �� � q-� ������
		void HermiteInterpolationFromCell(cell& q);

		/// ������ �������� ���������
		void CalcTotalVelo();

		/// ������ ���������� � 16 ��������� ����� � ������� ��������� ��������
		numvector<double, 16> CalcGam16(const cell& q);

		/// ����������
		~FastFourier() {};
	};

}//namespace FFT

#endif