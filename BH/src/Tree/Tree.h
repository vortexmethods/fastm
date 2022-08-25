/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.1    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/08/24     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2022 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Tree.h                                                           |
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
\brief Заголовок класса Tree
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.1
\date 24 августа 2022 г.
*/

#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <memory>

#include "Params.h"
#include "PointsCopy.h"

namespace BH
{

	extern long long op;


	struct TParticleCode
	{
		unsigned int key;
		int originNumber;
		Point2D r;
		double g;
	};


	//Структура ячейки мортоновского дерева 
	struct treeCellT
	{
		int parent;
		numvector<int, 2> child; //потомки: 0 - левый, 1 - правый
		numvector<int, 2> range; //диапазон частиц: 0 - первая частица в ячейке, 1 - последняя частица в ячейке
		bool particle; //признак того, что данная вершина - частица, а не кластер 
		               //(нужно не для обхода дерева, а для "адресации" частиц и кластеров)

		bool leaf;     //признак того, что глубже дерево обходить не нужно

		//длины префиксов и сами префиксы - в явном виде не нужны
		//int prefixLength;
		//unsigned int prefix;

		Point2D size;   //размер ячейки по x и по y (для листа (0;0) )
		Point2D centre; //центр ячейки

		int level;      //уровень текущей ячейки в дереве

		/*
		// Монопольный момент ячейки 
		double mon;
		// Дипольный, квадрупольный, октупольный и гексадекапольный моменты ячейки 
		Point2D dip, qua, oct, hex;
		*/
		
		//Мультипольные моменты ячейки 
		numvector<Point2D, order + 1> mom;

		//Коэффициенты для вычисления скоростей
		numvector<Point2D, order> E;

		//Вектор указателей на ячейки в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для нижних уровней
		std::vector<int> closeCells;
	};


	//Мортоновское дерево
	class MortonTree
	{
	private:
		//Ссылка на список обобщенных вихрей, по которым строится дерево
		std::vector<PointsCopy>& pointsCopy;
		
		//Список мортоновских кодов, вычисляемых по положениям частиц
		std::vector<TParticleCode> mortonCodes;

		//Внутренние временные переменные, нужны для сортировки
		std::vector<TParticleCode> mortonCodes_temp;
		std::vector<unsigned int> s;


		//Функция вычисления длины общей части префиксов двух (а значит - и диапазона) частиц
		int delta(int i, int j);

		//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
		int prefixLength(int cell);

		//Функция вычисления общего префикса двух частиц
		std::pair<unsigned int, int> prefix(int cell);
		
		//Функция вычисления геометрических параметров внутренней ячейки дерева
		void setCellGeometry(int cell);

		//Масштабный фактор дерева
		double iQuadSideVar;

		//Нижний левый угол дерева
		Point2D lowLeft;

		//Функция построения i-й внутренней ячейки дерева
		void buildInternalTreeCell(int i);
		
		//Основная функция вычисления вихревого влияния
		void influenceComputation(std::vector<Point2D>& result);
		
		//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
		static unsigned int expandBits(unsigned int v);

		//Мортоновский код для пары из чисел типа double
		//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
		static unsigned int morton2D(const Point2D& r);

		//Поиск габаритного прямоугольника системы точек
		//работает под OpenMP 2.0
		static std::pair<Point2D, Point2D> find_enclosing_rectangle_old(const std::vector<PointsCopy>& points);
		
		//Поиск габаритного прямоугольника системы точек
		//работает под OpenMP 4.0
		//static std::pair<Point2D, Point2D> find_enclosing_rectangle(const std::vector<PointsCopy>& points);

		//Поразрядная сортировка мортоновских кодов
		void RSort_Node3(TParticleCode* m, TParticleCode* m_temp, unsigned int n);
		//Вспомогательная функция к предыдущей
		inline void RSort_step3(TParticleCode* source, TParticleCode* dest, unsigned int n, unsigned int* offset, unsigned char sortable_bit);
		
		//Параллельная версия поразрядной сортировки мортоновских кодов
		inline void RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s);

	public:
		numvector<double, order + 1> iFact;
		numvector<double, (order + 1) * (order + 1)> binomCft;
		double& getBinom(int n, int k) {
			return binomCft[n * (order + 1) + k];
		}
		
		
		//Конструктор
		MortonTree(std::vector<PointsCopy>& points);		
		
		//Вектор номеров ячеек нижнего уровня 
		std::vector<int> mortonLowCells;

		//Дерево
		//Размер дерева равен сумме числа внутренних вершин (число частиц без единицы) и числа частиц
		//резервируется на 1 штуку больше для удобства работы
		//номера [0 ... pos.size()-2] отвечают внутренним узлам,
		//номер pos.size()-1 не используется 
		//номера [pos.size() ... 2*pos.size()-1] отвечают частицам
		std::vector<treeCellT> mortonTree;

		//Построение корня дерева и задание его общих параметров
		void makeRootMortonTree();

		//Построение внутренних ячеек дерева
		void buildMortonInternalTree();

		//Построение верхушек дерева --- отдельных частиц
		void buildMortonParticlesTree();

		//Два варианта заполнения списка нижних вершин: 
		//первый рекурсивный, быстро работает в последовательном варианте
		//второй нерекурсивный, хорошо распараллеливается, но неэффективен в последовательном варианте
		void fillMortonLowCells(int cell = 0);
		void fillMortonLowCellsA();

		//Расчет мультипольных моментов всех ячеек дерева
		void calculateMortonTreeParams(int cell = 0);
		
		//Расчет коэффицентов локального разложения для нижних ячеек 
		//calcCloseTrees --- признак заполнения списка ближных ячеек (нижнего уровня) к ячейке lowCell
		void calcLocalCoeffToLowLevel(int lowCell, int fromWho, bool calcCloseTrees);		
		
		//Расчет влияния от ближних ячеек по формуле Био - Савара
		void calcVeloBiotSavart(int lowCell);

		//Расчет влияния от дальней зоны при помощи Тейлоровского разложения
		void calcVeloTaylorExpansion(int lowCell);
	};


	//умножение комплексных чисел
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[1] * b[0] + a[0] * b[1] });
	}

	// умножение a на комплексно сопряженноe к b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
#ifdef calcOp
		op += 4;
#endif
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

	//возведение в квадрат
	template <typename T>
	inline T sqr(T x)
	{
		return x * x;
	}


	//Знак числа, написана оптимизированная версия, 
	//которая работает самым быстрым возможным образом, т.к. не содержит ни одного оператора условного перехода
	template <typename T>
	inline int sign(T val)
	{
		return (T(0) < val) - (val < T(0));
	}
	

	//Округление "в потолок" результата деления x на y
	//(годится для любых натуральных, но работает довольно медленно, здесь не нужна)
	//unsigned int ceil(unsigned int x, unsigned int y)
	//{
	//    return x / y + (x % y != 0);
	//}

	//Округление "в потолок" результата деления на степень двойки
	inline int ceilpow2(unsigned int x, unsigned int p) //  =ceil(x / 2^p)
	{
		return (x >> p) + !!(x & ((1 << p) - 1));
	}

	//Округление "в потолок" результата деления пополам
	inline int ceilhalf(unsigned int x) //  =ceil(x / 2), т.е. предыдущая функция при p=1
	{
		return (x >> 1) + (x & 1);
	}

}//namespace BH


#endif
