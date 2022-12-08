/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.3    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2022/12/08     |
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
\version 1.3
\date 08 декабря 2022 г.
*/

#ifndef TREE_H_
#define TREE_H_


#include <memory>

#include "defs.h"

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
		Point2D center; //центр ячейки

		int level;      //уровень текущей ячейки в дереве
		
		//Мультипольные моменты ячейки 
		std::vector<Point2D> mom;

		//Коэффициенты для вычисления скоростей
		std::vector<Point2D> E;

		//Вектор указателей на ячейки в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для нижних уровней
		std::vector<int> closeCells;
	};


	//Мортоновское дерево
	class MortonTree
	{
	private:
		//Ссылка на параметры
		const params& prm;

		//Ссылка на список обобщенных вихрей, по которым строится дерево
		std::vector<PointsCopy>& pointsCopy;

		//Список мортоновских кодов, вычисляемых по положениям частиц
		std::vector<TParticleCode> mortonCodes;

		//Внутренние временные переменные, нужны для сортировки
		std::vector<TParticleCode> mortonCodes_temp;
		std::vector<unsigned int> s;

		//Масштабный фактор дерева
		double iQuadSideVar;

		//Нижний левый угол дерева
		Point2D lowLeft;

		//Функция вычисления длины общей части префиксов двух (а значит - и диапазона) частиц
		int Delta(int i, int j) const;

		//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
		int PrefixLength(int cell) const;

		//Функция вычисления общего префикса двух частиц
		std::pair<unsigned int, int> Prefix(int cell) const;
		
		//Функция вычисления геометрических параметров внутренней ячейки дерева
		void SetCellGeometry(int cell);

		//Функция построения i-й внутренней ячейки дерева
		void BuildInternalTreeCell(int i);
		
	
		//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
		static unsigned int ExpandBits(unsigned int v);

		//Мортоновский код для пары из чисел типа double
		//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
		static unsigned int Morton2D(const Point2D& r);

		//Поиск габаритного прямоугольника системы точек
		//работает под OpenMP 2.0
		static std::pair<Point2D, Point2D> FindEnclosingRectangleOld(const std::vector<PointsCopy>& points);
		
		//Поиск габаритного прямоугольника системы точек
		//работает под OpenMP 4.0
		//static std::pair<Point2D, Point2D> FindEnclosingRectangle(const std::vector<PointsCopy>& points);

		//Поразрядная сортировка мортоновских кодов
		void RSort_Node3(TParticleCode* m, TParticleCode* m_temp, unsigned int n);

		//Вспомогательная функция к предыдущей
		inline void RSort_step3(TParticleCode* source, TParticleCode* dest, unsigned int n, unsigned int* offset, unsigned char sortable_bit);
		
		//Параллельная версия поразрядной сортировки мортоновских кодов
		inline void RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s);

	public:
		std::vector<double> iFact;
		std::vector<double> binomCft;
		
		void getStatistics(int& numParticles, int& treshold, std::pair<int, int>& partLevel, std::pair<int, int>& leafsLevel, int& lowLevelCells) const;
		
		double& getBinom(int n, int k) {
			return binomCft[n * (prm.order + 1) + k];
		}
		
		//std::vector<double> tempBuffer;
		
		//Признак того, что дерево из панелей (иначе - из точек)
		bool pans;

		//Максимальная глубина дерева при его обходе
		const int maxTreeLevel;

		//Конструктор
		MortonTree(const params& prm_, int maxTreeLevel_, std::vector<PointsCopy>& points, bool ifpans);
		
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
		void MakeRootMortonTree();

		//Построение внутренних ячеек дерева
		void BuildMortonInternalTree();

		//Построение верхушек дерева --- отдельных частиц
		void BuildMortonParticlesTree();

		//Два варианта заполнения списка нижних вершин: 
		//первый рекурсивный, быстро работает в последовательном варианте
		//второй нерекурсивный, хорошо распараллеливается, но неэффективен в последовательном варианте
		void FillMortonLowCells(int cell = 0);
		void FillMortonLowCellsA();

		///Расчет мультипольных моментов всех ячеек дерева
		// \param[in] cell --- ячейка для вычисления мультипольных моментов 
		void CalculateMortonTreeParams(int cell, int omplevel);
		
		//Расчет коэффицентов локального разложения для нижних ячеек 
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой будет рассчитываться влияние
		// \param[in] treeInf --- ссылка на влияющее дерево
		// \param[in] fromWho --- индекс влияющей ячейки
		// \param[in] calcCloseTrees --- признак заполнения списка ближных ячеек (нижнего уровня) к ячейке lowCell
		void CalcLocalCoeffToLowLevel(int lowCell, std::unique_ptr<MortonTree>& treeInf, int fromWho, bool calcCloseTrees);
		
		///Расчет влияния от ближних ячеек по формуле Био - Савара
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние//\param[in] treeInf --- влияющее дерево (чтобы брать ближние ячейки по индексу из нужного массива)
		// \param[in] treeInf --- влияющее дерево (чтобы брать ближние ячейки по индексу из нужного массива)
		// \param[in] type --- тип объектов во влияющем дереве: 
		// type = 0 для вихрей
		// type = 1 для панелей 
		void CalcVeloBiotSavart(int lowCell, std::unique_ptr<MortonTree>& treeInf);
		//void CalcVeloBiotSavartType0(int lowCell, std::unique_ptr<MortonTree>& treeInf);

		///Расчет интегрального влияния от распределения завихренности панелей ближних ячеек (для решения СЛАУ)
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void CalcInfluenceFromPanels(int lowCell);

		/// Обновление влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void UpdateInfluence(int lowCell);

		///Расчет точечного влияния от распределения завихренности панелей ближних ячеек (для решения СЛАУ)
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void CalcInfluenceFromPanelsToPoints(int lowCell);

		///Расчет влияния от следа на панели профиля
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		// \param[in] treeInf --- влияющее дерево (чтобы брать ближние ячейки по индексу из нужного массива)
		void CalcInfluenceFromPointsToPanel(int lowCell, std::unique_ptr<MortonTree>& treeInf);

		///Расчет влияния от дальней зоны при помощи Тейлоровского разложения
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void CalcVeloTaylorExpansion(int lowCell);

		///Добавление влияния дальней зоны для правой части
		// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void AddVelo(int lowCell);
	};
	   
}//namespace BH


#endif
