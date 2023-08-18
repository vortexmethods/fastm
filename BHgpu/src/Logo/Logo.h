/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.5    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/08/29     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: Logo.h                                                           |
| Info: Source code of BHgpu                                                  |
|                                                                             |
| This file is part of BHgpu.                                                 |
| BHcu is free software: you can redistribute it and/or modify it             |
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
| along with BHgpu.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Логотип
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#ifndef LOGO_H_
#define LOGO_H_

#include <fstream>
#include "Params.h"

namespace BHcu
{

	void PrintLogoToStream(std::ostream& str)
	{
		str <<
			"/*--------------------------------*- BHgpu -*----------------*---------------*\\" << '\n' << \
			"| #####   ##  ##                |                            | Version 1.5    |" << '\n' << \
			"| ##  ##  ##  ##   ####  ##  ## |  BHcu: Barnes-Hut method   | 2023/08/29     |" << '\n' << \
			"| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*" << '\n' << \
			"| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |" << '\n' << \
			"| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |" << '\n' << \
			"| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |" << '\n' << \
			"\\*---------------------------------------------------------------------------*/" << '\n';
	}


	void PrintConfiguration(int nbodies)
	{
		std::cout << std::endl;
		std::cout << "                     Configuration for task \"" << task << "\"" << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		std::cout << "Floating point type:           " << (sizeof(real) == 4 ? "float" : "double") << " (" << sizeof(real) << " bytes)" << std::endl;
		std::cout << "Number of bodies:              " << nbodies << " (eps = " << EPS << ")" << std::endl;
		std::cout << "Multipoles:                    " << "up to " << order << "-th order" << std::endl;
		std::cout << "Theta in proximity criterion:  " << THETA << std::endl;
		std::cout << "Compare/save:                  " << (compare ? "true" : "false") << "/" << (save ? "true" : "false") << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}



	void PrintStatistics(int run, int runs,
		const float* timing, const float* mintiming, const float* avtiming,
		double runtime, double minruntime, double avruntime)
	{
		if (run == 0)
		{
			std::cout << std::endl;
			std::cout << "                         Time statistics                                " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "  run/runs  total      BBK   TBK   CLK   SKK     FCK      ker.time      " << std::endl;
		}


		printf("# %3d/%-3d: %6.3lf s  (", run + 1, runs, runtime);

		printf(" %3.1f ", timing[1]);
		printf(" %4.1f ", timing[2]);
		printf(" %4.1f ", timing[3]);
		printf(" %4.1f ", timing[4]);
		printf(" %6.1f ", timing[5]);

		printf(") = %6.1f ms\n", timing[6]);		

		if (run == runs - 1)
		{
			std::cout << "------------------------------------------------------------------------" << std::endl;

			printf("     min : %6.3lf s  (", minruntime);

			printf(" %3.1f ", mintiming[1]);
			printf(" %4.1f ", mintiming[2]);
			printf(" %4.1f ", mintiming[3]);
			printf(" %4.1f ", mintiming[4]);
			printf(" %6.1f ", mintiming[5]);
			printf(") | %6.1f ms\n", mintiming[6]);


			printf("    aver : %6.3lf s  (", avruntime / runs);

			printf(" %3.1f ", avtiming[1] / runs);
			printf(" %4.1f ", avtiming[2] / runs);
			printf(" %4.1f ", avtiming[3] / runs);
			printf(" %4.1f ", avtiming[4] / runs);
			printf(" %6.1f ", avtiming[5] / runs);
			printf(") | %6.1f ms\n", avtiming[6] / runs);
		}
	}


	void PrintAccuracyHead()
	{
		std::cout << std::endl;
		std::cout << "                         Accuracy control                               " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	void PrintAccuracyError(real val)
	{
		std::cout << "Relative error:     " << val << std::endl;

		std::cout << "------------------------------------------------------------------------" << std::endl;
	}



	
}
#endif