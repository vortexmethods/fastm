/*--------------------------------*- FFT -*------------------*---------------*\
|    ######  ######  ######     |                            | Version 1.0    |
|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |
|    ####    ####      ##       |  for 2D vortex particles   *----------------*
|    ##      ##        ##       |  Open Source Code                           |
|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
*-----------------------------------------------------------------------------*
| File name: Logo.h                                                           |
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
\brief Логотип
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\version 1.0
\date 05 августа 2021 г.
*/

#ifndef LOGO_H_
#define LOGO_H_

#include <fstream>

namespace FFT
{

	void PrintLogoToStream(std::ostream& str)
	{
		str <<
			"/*--------------------------------*- FFT -*------------------*---------------*\\" << '\n' << \
			"|    ######  ######  ######     |                            | Version 1.0    |" << '\n' << \
			"|    ##      ##        ##       |  FFT: FFT-based method     | 2021/08/05     |" << '\n' << \
			"|    ####    ####      ##       |  for 2D vortex particles   *----------------*" << '\n' << \
			"|    ##      ##        ##       |  Open Source Code                           |" << '\n' << \
			"|    ##      ##        ##       |  https://www.github.com/vortexmethods/fastm |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |" << '\n' << \
			"\\*---------------------------------------------------------------------------*/" << '\n';
	}

	void PrintProperties()
	{
		printf("\n");
		printf("                             CPU Properties                             \n");
		printf("------------------------------------------------------------------------\n");
		printf("OpenMP max threads:                           %d                        \n", omp_get_max_threads());
		printf("Number of OpenMP threads for TreeBuilding:    %d                        \n", numTh1);
		printf("Number of OpenMP threads for VeloComputing:   %d                        \n", numTh3);
		printf("------------------------------------------------------------------------\n");
	}



	void PrintConfiguration(int nbodies)
	{
		std::cout << std::endl;
		std::cout << "                     Configuration for task \"" << task << "\"" << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		std::cout << "Number of bodies:               " << nbodies << " (eps = " << eps << ")" << std::endl;
		std::cout << "Mesh:                           " << NP << " x " << NP << std::endl;
		std::cout << "Size of close zone:             " << sizeOfCloseZone << " layers" << std::endl;
		std::cout << "Interpolation mode:             ";
		if (interpType == 0)
			std::cout << "bi-linear" << std::endl;
		if (interpType == 1)
			std::cout << "bi-hermitian" << std::endl;		
		std::cout << "Compare/save:                   " << (compare ? "true" : "false") << "/" << (save ? "true" : "false") << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}


	void PrintStatistics(int run, int runs,
		const double* timing, const double* mintiming, const double* avtiming,
		double runtime, double minruntime, double avruntime)
	{
		if (run == 0)
		{
			std::cout << std::endl;
			std::cout << "                         Time statistics                                " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "  run/runs  total        MSH     COR     FFT     VEL      ker.time      " << std::endl;
		}


		printf("# %3d/%-3d: %6.3lf s  (", run + 1, runs, runtime);

		printf(" %6.3f ", timing[1]);
		printf(" %6.3f ", timing[2]);
		printf(" %6.3f ", timing[3]);
		printf(" %6.3f ", timing[4]);
		printf(") = %6.3f s\n", timing[5]);

		if (run == runs - 1)
		{
			std::cout << "------------------------------------------------------------------------" << std::endl;

			printf("     min : %6.3lf s  (", minruntime);

			printf(" %6.3f ", mintiming[1]);
			printf(" %6.3f ", mintiming[2]);
			printf(" %6.3f ", mintiming[3]);
			printf(" %6.3f ", mintiming[4]);
			printf(") | %6.3f s\n", mintiming[5]);


			printf("    aver : %6.3lf s  (", avruntime / runs);

			printf(" %6.3f ", avtiming[1] / runs);
			printf(" %6.3f ", avtiming[2] / runs);
			printf(" %6.3f ", avtiming[3] / runs);
			printf(" %6.3f ", avtiming[4] / runs);
			printf(") | %6.3f s\n", avtiming[5] / runs);
			/*
			#ifdef calcOp
						std::cout << "------------------------------------------------------------------------" << std::endl;
						std::cout << "Operations:   " << op << " (parallelization is partially turned off)    " << std::endl;
			#endif
			*/
		}
	}






	void PrintAccuracyHead()
	{
		std::cout << std::endl;
		std::cout << "                         Accuracy control                               " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	void PrintAccuracyError(double val)
	{
		std::cout << "Relative error:     " << val << std::endl;

		std::cout << "------------------------------------------------------------------------" << std::endl;
	}


}//namespace FFT
#endif