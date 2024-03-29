#/*---------------------------------*- BH -*------------------*---------------*\
#|        #####   ##  ##         |                            | Version 1.3    |
#|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
#|        #####   ######         |  for 2D vortex particles   *----------------*
#|        ##  ##  ##  ##         |  Open Source Code                           |
#|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
#|                                                                             |
#| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
#*-----------------------------------------------------------------------------*
#| File name: CMakeLists.txt                                                   |
#| Info: Source code of BH                                                     |
#|                                                                             |
#| This file is part of BH.                                                    |
#| BH is free software: you can redistribute it and/or modify it               |
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
#| along with BH.  If not, see <http://www.gnu.org/licenses/>.                 |
#\*---------------------------------------------------------------------------*/

cmake_minimum_required(VERSION 3.0)
project(BH LANGUAGES CXX)
#set(CMAKE_CXX_STANDARD 17)



message (STATUS "C++ compiler is: " ${CMAKE_CXX_COMPILER})


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

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

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
find_package(OpenMP REQUIRED)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	
#message (STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")
message (STATUS "OpenMP version: ${OpenMP_CXX_VERSION_MAJOR}")
IF(OpenMP_CXX_VERSION_MAJOR LESS 3.0)
	add_definitions(-DOLD_OMP)
endif()

OPTION(NBODY "N-body problem, computation velo of vortex particles" OFF)
IF(NBODY)
    ADD_DEFINITIONS(-DCALCVORTEXVELO)
	message (STATUS "N-body problem, computation velo of vortex particles")
ENDIF(NBODY)

OPTION(BIE "BIE solution, vortex sheet intensity comptation" OFF)
IF(BIE)
    ADD_DEFINITIONS(-DCALCSHEET)
	message (STATUS "BIE solution, vortex sheet intensity comptation")
ENDIF(BIE)

OPTION(VP "Velo calculation in fluid domain" OFF)
IF(VP)
    ADD_DEFINITIONS(-DCALCVP)
	message (STATUS "Velo calculation in fluid domain")
ENDIF(VP)



OPTION(LINEAR "Piecewise-linear solution representation over panels" OFF)
IF(LINEAR)
    ADD_DEFINITIONS(-DlinScheme)
	message (STATUS "Piecewise-linear solution representation over panels")
ENDIF(LINEAR)

OPTION(AS "Piecewise-linear solution representation with asymptotics resolving over panels" OFF)
IF(AS)
    ADD_DEFINITIONS(-DlinScheme)
	ADD_DEFINITIONS(-DasympScheme)
	message (STATUS "Piecewise-linear solution representation with asymptotics resolving over panels")
ENDIF(AS)




OPTION(CALCOP "Count operations" OFF)
IF(CALCOP)
    ADD_DEFINITIONS(-DcalcOp)
	message (STATUS "Counting operations")
ENDIF(CALCOP)


ADD_DEFINITIONS(-D_SILENCE_CXX17_NEGATORS_DEPRECATION_WARNING)
ADD_DEFINITIONS(-D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})