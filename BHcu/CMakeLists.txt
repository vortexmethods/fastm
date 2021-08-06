#/*--------------------------------*- BHcu -*-----------------*---------------*\
#| #####   ##  ##                |                            | Version 1.0    |
#| ##  ##  ##  ##   ####  ##  ## |  BHcu: Barnes-Hut method   | 2021/08/05     |
#| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
#| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
#| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
#|                                                                             |
#| Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina                   |
#| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
#*-----------------------------------------------------------------------------*
#| File name: CMakeLists.txt                                                   |
#| Info: Source code of BHcu                                                   |
#|                                                                             |
#| This file is part of BHcu.                                                  |
#| BHcu is free software: you can redistribute it and/or modify it             |
#| under the terms of the GNU General Public License as published by           |
#| the Free Software Foundation, either version 3 of the License, or           |
#| (at your option) any later version.                                         |
#|                                                                             |
#| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
#| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
#| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
#| for more details.                                                           |
#|                                                                             |
#| You should have received a copy of the GNU General Public License           |
#| along with BHcu.  If not, see <http://www.gnu.org/licenses/>.               |
#\*---------------------------------------------------------------------------*/

cmake_minimum_required(VERSION 2.8)
project(BHcu LANGUAGES CXX CUDA)

enable_language(CUDA)
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_62,code=sm_62 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_72,code=sm_72 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75")

set(CMAKE_CONFIGURATION_TYPE "Debug;Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(SRC_LIST src/main/main.cpp)

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cuh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cuh
  )

add_executable(${PROJECT_NAME} ${ALLFILES})

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${dirlist})
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${dirlist})

#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
endif()

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})