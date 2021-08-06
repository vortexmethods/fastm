fastm: Fast algorithms for 2D vortex particle method
====================================================

<p align="center"><img src="https://www.mdpi.com/entropy/entropy-23-00118/article_deploy/html/images/entropy-23-00118-g019-550.jpg"></p>

![Repo Size](https://img.shields.io/github/repo-size/vortexmethods/fastm.svg)
![License](https://img.shields.io/github/license/vortexmethods/fastm.svg)

Copyright (C) 2020-2021 Ilia Marchevsky, Evgeniya Ryatina

Программная реализация (с открытым исходным кодом) быстрых алгоритмов расчета скоростей вихревых частиц для двумерных методов вычислительной гидродинамики.

Версия 1.0 от 05 августа 2021 г.

ЛИЦЕНЗИЯ
--------

Программа распространяется на условиях свободной лицензии [GNU GPLv3](https://www.gnu.org/licenses/gpl.txt)
   
   
ИСПОЛЬЗУЕМЫЕ МЕТОДЫ
-------------------
   
В репозитории представлены программные реализации следующих методов:

* ***BH*** - [метод Барнса-Хата, BH](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation), оригинальная версия которого была описана [Josh Barnes и Piet Hut в 1986 г.](https://www.nature.com/articles/324446a0), представленная здесь модификация разработана авторами; 

* ***BHcu*** - [CUDA-реализация](http://cs.txstate.edu/~mb92/papers/gcg11.pdf) модификации метода Барнса-Хата, разработанная авторами в результате переработки и развития кода, первоначально разработанного [Martin Burtscher](https://userweb.cs.txstate.edu/~burtscher/) и доступного под свободной лицензией *Modified BSD License* на [GitHub](https://github.com/IntelligentSoftwareSystems/Galois/tree/master/lonestar/scientific/gpu/barneshut); 

* ***FMM*** - [быстрый метод мультиполей](https://ru.wikipedia.org/wiki/%D0%91%D1%8B%D1%81%D1%82%D1%80%D1%8B%D0%B9_%D0%BC%D0%B5%D1%82%D0%BE%D0%B4_%D0%BC%D1%83%D0%BB%D1%8C%D1%82%D0%B8%D0%BF%D0%BE%D0%BB%D0%B5%D0%B9), оригинальная версия которого была описана [Leslie Greengard и Владимиром Рохлиным (мл.) в 1987 г.](https://www.sciencedirect.com/science/article/pii/0021999187901409), представленная здесь реализация разработана авторами совместно с Д.О. Попудняк;  

* ***FFT*** - метод на основе решения задачи Пуассона для функции тока методом быстрого преобразования Фурье, оригинальная версия которого была описана [Guido Morgenthal и Jens Honore Walther в 2007 г.](https://www.sciencedirect.com/science/article/abs/pii/S004579490700034X); представленная здесь реализация со схемами интерполяции повышенной точности разработана авторами; 


СТРУКТУРА РЕПОЗИТОРИЯ
---------------------

* *BH* - исходные коды метода Барнса-Хата в авторской модификации

* *BHcu* - исходные коды CUDA-реализации модифицированного метода Барнса-Хата

* *FMM* - исходные коды быстрого метода мультиполей

* *FFT* - исходные коды метода на основе быстрого преобразования Фурье

* *include* - общие библиотеки (eigen и вспомогательные типы данных)

* *test* - примеры исходных данных --- распределения различного количества вихревых частиц в единичном квадрате (**файлы с большим количеством вихрей заархивированы для экономии места!**)

* *res* - папка, в которую производится сохранение результатов; в нее помещены результаты прямого решения тестовых задач методом Био-Савара (**файлы для большого количества вихрей заархивированы для экономии места!**)


	 
УСТАНОВКА
---------

Для запуска программ *BH*, *BHcu*, *FMM*, *FFT* на комьютере необходимо загрузить исходные коды программы.
Если на Вашем компьютере установлен "Git", достаточно исполнить команду 

      git clone https://github.com/vortexmethods/fastm.git fastm

по результатам работы которой в текущей папке будет создана подпапка fastm и в нее будут загружены все файлы из репозитория.	  
	  
Подготовка к компияции исходных кодов предполагает создание папок "build" в каталогах *BH*, *BHcu*, *FMM*, *FFT* с загруженными исходными кодами, переход в эти папки и выполнение команды 

      cmake ..
	  
При необходимости следует указать необходимые ключи для настройки используемых компиляторов, указания опций компиляции и т.п., возможно, потребуется также некоторая модификация файла "CMakeLists.txt", содержащего параметры настройки CMake.

В частности, для подготовки исходных кодов для их последующей компиляции в Windows средствами MS Visual Studio следует, в зависимости от версии, использовать одну из следующих команд (опция Win64 обязательна для использования возможности проведения вычислений на графических картах Nvidia CUDA, для Visual Studio 2019 она включена по умолчанию для 64-битных систем)
      
      cmake -G"Visual Studio 15 2017 Win64" ..
      cmake -G"Visual Studio 16 2019" ..
	  

В случае использования компилятора, отличного от используемого по умолчанию в Windows (это, как правило, встроенный в MS Visual Studio компилятор MVSC), например, компилятора Intel, при подготовке исходных кодов к компиляции необходимо указать, в зависимости от версии, ключ (отметим, что Intel C++ Compiler 19 интегрируется в Visual Studio 2019 лишь начиная с версии Upd.4)
      
	  -T"Intel C++ Compiler <ver>"

При работе в Linux альтернативный компилятор с C++ (к примеру, icpc для компилятора Intel вместо используемого в большинстве случаев по умолчанию компилятора g++) требуется исполнить команду

      CXX=icpc cmake ..

Дальнейшая компиляция кода зависит от используемой операционной системы. В Windows, как правило, при помощи CMake будет создан проект для его поледующего открытия и компиляции средствами MS Visual Studio (см. выше), в Linux достаточно исполненияиз созданной папки команды

      make

	 
РАБОТА С CUDA
-------------
	 
Запуск программы *BHcu* возможен только при наличии в системе установленного [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit). В ходе работы cmake программа будет автоматически сконфигурирована для счета с использованием возможностей графических карт, работающих по технологии NVidia CUDA.
 
При выполнении расчетов на узле, на котором установлено несколько видеокарт, по умолчанию все расчеты будут производиться на устройстве с индексом #0. Чтобы этого избежать, нужно открыть файл исходного кода BHcu/src/Params/Params.h и указать соответствующий индекс в качестве значения переменной dev.
	  
	  
ЗАПУСК
------	  
	  	  
Запуск программ осуществляется строго из папок build, что необходимо для коректного считывания исходных файлов (из папки test) и сохранения результатов (в папку res)


НЕОБХОДИМОЕ ПО
--------------

Для компиляции требуется наличие установленных: 

* системы автоматизации сборки программного обеспечения из исходного кода [cmake](https://cmake.org/),
* компилятора с языка C++, поддерживающего технологию OpenMP и стандарт С++11,
* библиотеки [Eigen](http://eigen.tuxfamily.org) (не обязательно, исходные коды находятся в папке "include"),
* при наличии в системе графического процессора Nvidia, поддерживающего технологию [CUDA](https://ru.wikipedia.org/wiki/CUDA), для его использования необходим [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit).


КРАТКАЯ ИСТОРИЯ ВЕРСИЙ
----------------------

* Версия 1.0 представлена 09 августа 2021 г. в рамках программы ["Вычислительные технологии, многомерный анализ данных и моделирование"](https://sochisirius.ru/obuchenie/graduates/smena901/4333) (Образовательный центр "Сириус", Сочи, 2-22 августа 2021 г.)


ВОПРОСЫ, ПРЕДЛОЖЕНИЯ И ЗАМЕЧАНИЯ
--------------------------------

На [странице Issues](https://github.com/vortexmethods/fastm/issues) мы будем рады ответить на Ваши вопросы, с благодарностью выслушаем предложения и замечания по коду.

Мы открыты для любого конструктивного взаимодействия!


---
С глубоким уважением,
разработчики

<table width="500" border="0" cellpadding="5">

<tr>

<td align="center" valign="center" width="20%">
<img src="https://raw.githubusercontent.com/vortexmethods/VM2D/master/docs/_static/authors/Marchevsky.jpg" alt="Марчевский И.К."/>
<br />
Марчевский Илья Константинович
</td>

<td align="center" valign="center" width="20%">
<img src="https://raw.githubusercontent.com/vortexmethods/VM2D/master/docs/_static/authors/Ryatina.jpg" alt="Рятина Е.П."/>
<br />
Рятина Евгения Павловна
</td>

</tr>

</table>






