/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
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
\brief ��������� �������� ������
\author ���������� ���� ��������������
\author ������ ������� ��������
\version 1.0
\date 05 ������� 2021 �.
*/


#ifndef PARAMS_H_
#define PARAMS_H_

namespace FFT
{

#include <iostream>
	// ������ ��������� ��������
	static const double eps = 1e-3;

	// ��� ����� � �������
	static const std::string nameFile = "../../test/test100000.txt";
	// �������� ������ ��� ����� � �����������
	static const std::string task = "100k";

	// ����� ����� �����
	static const int NP = 128;
	// ���������� ����������� ����� � ������� ����
	static const int sizeOfCloseZone = 3;
	// ��� ������������: 0 - ����������, 1 - ��������
	static const int interpType = 1;

	// ���������� OMP ����� ��� 1 � 3 �������� �������������� (0 - ����������� ���������)
	static const int numTh1 = 4;
	static const int numTh3 = 4;

	//��������� � ������ �������
	static const bool compare = true;

	//����� ��������� �������� ����
	static const int runs = 10;


	//��� �������, ��������� ���� �� ����� ������

	// ������ ��������� �������� (� ��������)
	static const double eps2 = eps * eps;

	//������ ��������� �� ������� ������ �� ����� (���� ��� ���, �� ����� �������� � ��������)
	static const bool BSfromFile = true;
	//���������� ��������� �� �������� ������ � ����
	static const bool save = false;

}//namespace FFT

#endif