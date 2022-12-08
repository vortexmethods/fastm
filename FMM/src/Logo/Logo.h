/*--------------------------------*- FMM -*------------------*---------------*\
|   ######  ##   ##  ##   ##    |                            | Version 1.3    |
|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |
|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*
|   ##      ##   ##  ##   ##    |  Open Source Code                           |
|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |
*-----------------------------------------------------------------------------*
| File name: Logo.h                                                           |
| Info: Source code of FMM                                                    |
|                                                                             |
| This file is part of FMM.                                                   |
| FMM is free software: you can redistribute it and/or modify it              |
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
| along with FMM.  If not, see <http://www.gnu.org/licenses/>.                |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Логотип
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Попудняк Дарья Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

#ifndef LOGO_H_
#define LOGO_H_

#include <fstream>

namespace FMM
{

	void PrintLogoToStream(std::ostream& str)
	{
		str <<
			"/*--------------------------------*- FMM -*------------------*---------------*\\" << '\n' << \
			"|   ######  ##   ##  ##   ##    |                            | Version 1.3    |" << '\n' << \
			"|   ##      ### ###  ### ###    |  FMM: Multipole method     | 2022/12/08     |" << '\n' << \
			"|   ####    ## # ##  ## # ##    |  for 2D vortex particles   *----------------*" << '\n' << \
			"|   ##      ##   ##  ##   ##    |  Open Source Code                           |" << '\n' << \
			"|   ##      ##   ##  ##   ##    |  https://www.github.com/vortexmethods/fastm |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2020-2022 Ilia Marchevsky, Evgeniya Ryatina, Daria Popudnyak  |" << '\n' << \
			"\\*---------------------------------------------------------------------------*/" << '\n';
	}


	void PrintProperties()
	{
		printf("\n");
		printf("                             CPU Properties                             \n");
		printf("------------------------------------------------------------------------\n");
		printf("OpenMP max threads:                 %d\n", omp_get_max_threads());		
		printf("------------------------------------------------------------------------\n");
	}

	void PrintConfiguration(int nbodies)
	{
		std::cout << std::endl;
		std::cout << "                     Configuration for task \"" << task << "\"" << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		std::cout << "Number of bodies:               " << nbodies << " (eps = " << eps << ")" << std::endl;
		std::cout << "Multipoles:                     " << nt << std::endl;
		std::cout << "Maximal tree level:             " << maxLevel << std::endl;
		std::cout << "Counting operations:            " <<
		#ifdef calcOp
					"true (OpenMP is partially turned off)" << std::endl;
		#else
					"false" << std::endl;
		#endif		
		std::cout << "Compare/save:                   " << (compare ? "true" : "false") << "/" << (save ? "true" : "false") << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}



	void PrintStatistics(int run, int runs,
		const double* timing, const double* mintiming, const double* avtiming,
		double runtime, double minruntime, double avruntime,
		long long op)
	{
		if (run == 0)
		{
			std::cout << std::endl;
			std::cout << "                         Time statistics                                " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "  run/runs  total       TREE    UPW    M2L    L2L    LEAF    ker.time   " << std::endl;
		}


		printf("# %3d/%-3d: %6.3lf s  (", run + 1, runs, runtime);

		printf(" %6.3f", timing[1]);
		printf(" %6.3f", timing[2]);
		printf(" %6.3f", timing[3]);
		printf(" %6.3f", timing[4]);
		printf(" %6.3f", timing[5]);
		printf(") = %6.3f s\n", timing[6]);

		if (run == runs - 1)
		{
			std::cout << "------------------------------------------------------------------------" << std::endl;

			printf("     min : %6.3lf s  (", minruntime);

			printf(" %6.3f", mintiming[1]);
			printf(" %6.3f", mintiming[2]);
			printf(" %6.3f", mintiming[3]);
			printf(" %6.3f", mintiming[4]);
			printf(" %6.3f", mintiming[5]);
			printf(") | %6.3f s\n", mintiming[6]);


			printf("    aver : %6.3lf s  (", avruntime / runs);

			printf(" %6.3f", avtiming[1] / runs);
			printf(" %6.3f", avtiming[2] / runs);
			printf(" %6.3f", avtiming[3] / runs);
			printf(" %6.3f", avtiming[4] / runs);
			printf(" %6.3f", avtiming[5] / runs);
			printf(") | %6.3f s\n", avtiming[6] / runs);
			#ifdef calcOp
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "Operations:   " << op << " operations                                   " << std::endl;
			#endif
			
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

}//namespace FMM
#endif