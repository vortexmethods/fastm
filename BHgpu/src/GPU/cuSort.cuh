/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.4    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/03/28     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: cuKernels.cu                                                     |
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
\brief Заголовочный файл для функций сортировки на GPU
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.4
\date 28 марта 2023 г.
*/

#ifndef CUSORT_H_
#define CUSORT_H_

#include <cuda.h>

namespace BHcu
{

    void RadixSortFromCUB(        
        const int* MmortonCodesKeyUnsortd, int* MmortonCodesKeyd, \
        const int* MmortonCodesIdxUnsortd, int* MmortonCodesIdxd,
        int num_items, int begin_bit, int end_bit);

    //(d_temp_storage, temp_storage_bytes, \
    //    MmortonCodesKeyUnsortd, MmortonCodesKeyd, \
    //    MmortonCodesIdxUnsortd, MmortonCodesIdxd, \
    //      nbodiesd, 0, 2 * codeLength);

}

#endif