/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Logo.h                                                           |
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
\brief Логотип
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.3
\date 08 декабря 2022 г.
*/

#ifndef LOGO_H_
#define LOGO_H_

#include <fstream>
#include <string>

#include "Params.h"

namespace BH
{
	/// Печать логотипа программы
	void PrintLogoToStream(std::ostream& str)
	{
		str <<
			"/*---------------------------------*- BH -*------------------*---------------*\\" << '\n' << \
			"|        #####   ##  ##         |                            | Version 1.3    |" << '\n' << \
			"|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |" << '\n' << \
			"|        #####   ######         |  for 2D vortex particles   *----------------*" << '\n' << \
			"|        ##  ##  ##  ##         |  Open Source Code                           |" << '\n' << \
			"|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |" << '\n' << \
			"\\*---------------------------------------------------------------------------*/" << '\n';
	}

	/// Печать числа используемых ядер при расчетах с помощью OpenMP и глубины порождения нитей при обходе дерева
	void PrintProperties(const params& prm)
	{
		printf("\n");
		printf("                             CPU Properties                             \n");
		printf("------------------------------------------------------------------------\n");
		printf("OpenMP max threads:                 %d\n", omp_get_max_threads());
		printf("Nested OpenMP up to tree level:     %d\n", prm.maxLevelOmp + 1);
		printf("------------------------------------------------------------------------\n");
	}


	/// \brief Печать текущей конфигурации
	///
	/// \param[in] prm константная ссылка на параметры, считанные из файла
	/// \param[in] nVrt число вихревых частиц в загруженном вихревом следе
	/// \param[in] nPnl число панелей на профиле
	/// \param[in] nVP число точек вычисления скорости в области течения
	void PrintConfiguration(const params& prm, int nVrt, int nPnl, int nVP)
	{
#ifdef CALCVORTEXVELO
		std::string label = "NBODY: ";
#else
#ifdef CALCSHEET
		std::string label = "  BIE: ";
#else
		std::string label = "   VP: ";
#endif
#endif
			
		auto shortName = [](const std::string& fileName)
		{
			char ch = '/';

			size_t index = fileName.rfind(ch);
			std::string shortAirfoilFile;
			if (index != std::string::npos)
				return fileName.substr(index + 1);
			else
				return fileName;
		};

		std::cout << std::endl;
		std::cout << "       " << label << " Configuration for task \"" << prm.task.c_str() << "\"" << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		std::cout << "Number of particles in wake:    " << nVrt << " (" << shortName(prm.vortexFile) << ")" << std::endl;
		std::cout << "Wake shift:                     " << prm.wakeShift << std::endl;
#if defined CALCSHEET || defined CALCVP
		std::cout << "Number of panels at airfoil:    " << nPnl << " (" << shortName(prm.airfoilFile) << ")" << std::endl;
#endif

#if defined CALCVP
		std::cout << "Number of VP points:            " << nVP << " (" << shortName(prm.vpFile) << ")" << std::endl;
#endif
		
		std::cout << "Epsilon:                        " << prm.eps << std::endl;
		std::cout << "Multipoles:                     " << "up to " << prm.order << "-th moment" << std::endl;
		std::cout << "Theta in proximity criterion:   " << prm.theta << std::endl;
		std::cout << "Counting operations:            " <<
#ifdef calcOp 
			"true (OpenMP is turned off; 'runs' = 1)" << std::endl;
#else
			"false" << std::endl;
#endif
		std::cout << "Compare/save:                   " << (prm.compare ? "true" : "false") << "/" << (prm.save ? "true" : "false") << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}


	/// \brief Печать информации о построенных деревьях
	///
	/// \param[in] treeVrt константная ссылка на дерево из вихревых частиц
	/// \param[in] treePan константная ссылка на дерево из центров панелей
	/// \param[in] treeVP константная ссылка на дерево из точек вычисления скоростей в области течения
	void PrintTreesInfo(const std::unique_ptr<MortonTree>& treeVrt, const std::unique_ptr<MortonTree>& treePan, const std::unique_ptr<MortonTree>& treeVP)
	{
		std::cout << std::endl;
		std::cout << "                         Trees statistics                               " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		std::cout << "Tree      N particles   Treshold   Part. depth   Leafs depth   N leafs  " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		
		//Печать статистики по одному дереву	
		auto printStat = [](const std::unique_ptr<MortonTree>& tree, const std::string& label)
		{
			if (tree.get() != nullptr)
			{
				int treshold, numParticles, cowLevelCells;
				std::pair<int, int> partLevel, leafsLevel;
				tree->getStatistics(numParticles, treshold, partLevel, leafsLevel, cowLevelCells);
				printf("%s  %7d        %2d       %2d ... %2d     %2d ... %2d   %7d\n", label.c_str(), numParticles, treshold, partLevel.first, partLevel.second, leafsLevel.first, leafsLevel.second, cowLevelCells);
			}
		};
		
		printStat(treeVrt, "Vortices: ");
		printStat(treePan, "Panels:   ");
		printStat(treeVP,  "VP points:");		
	}//PrintTreesInfo
		
	
	///Печать общей временнОй статистики
	void PrintStatistics(int run, int runs,
		const double* timing, const double* mintiming, const double* avtiming,
		double runtime, double minruntime, double avruntime, int niter)
	{
		if (run == 0)
		{
			std::cout << std::endl;
			std::cout << "                         Time statistics                                " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "  run/runs  total        TBK     SKK     FCK    | ker.time              " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
		}


		printf("# %3d/%-3d: %6.3lf s  (", run + 1, runs, runtime);

		printf(" %6.3f ", timing[1]);
		printf(" %6.3f ", timing[2]);
		printf(" %6.3f ", timing[3]);
		printf(") | %6.3f s", timing[4]);
		if (niter > 0)
			printf("   (%d iters)", niter);
		printf("\n");
		
		if (run == runs - 1)
		{
			std::cout << "------------------------------------------------------------------------" << std::endl;

			printf("     min : %6.3lf s  (", minruntime);

			printf(" %6.3f ", mintiming[1]);
			printf(" %6.3f ", mintiming[2]);
			printf(" %6.3f ", mintiming[3]);
			printf(") | %6.3f s\n", mintiming[4]);

			printf("    aver : %6.3lf s  (", avruntime / runs);

			printf(" %6.3f ", avtiming[1] / runs);
			printf(" %6.3f ", avtiming[2] / runs);
			printf(" %6.3f ", avtiming[3] / runs);
			printf(") | %6.3f s\n", avtiming[4] / runs);
			std::cout << "------------------------------------------------------------------------" << std::endl;
		}
	}

	/// Печать количества операций
	void PrintOps(long long nops)
	{
#ifdef calcOp
		if (nops > 0)
		{
			std::cout << "Operations:   " << nops << "  ( " << (double)(nops) << " )" << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
		}
#endif
	}

	/// Печать заголовка контроля точности
	void PrintAccuracyHead()
	{
		std::cout << std::endl;
		std::cout << "                         Accuracy control                               " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	/// Печать величины ошибки
	void PrintAccuracyError(double val)
	{
		std::cout << "Relative error:     " << val << std::endl;

		std::cout << "------------------------------------------------------------------------" << std::endl;
	}


}
#endif