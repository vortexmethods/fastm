#include "WolframLibrary.h"
#include "WolframCompileLibrary.h"

#include "numvector.h"
#include "source.h"

#include <vector>

using namespace std;

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion(){
	return WolframLibraryVersion;
}

EXTERN_C DLLEXPORT int WolframLibrary_initialize(WolframLibraryData \
	libData) {
	return 0;
}

EXTERN_C DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData \
	libData) {
	return;
}

//Возвращает матрицу скосов 
EXTERN_C DLLEXPORT int Matrc(WolframLibraryData libData,
	mint Argc, MArgument *Args, MArgument Res)
{
	double* pts = (double*)(MArgument_getMTensor(Args[0])->data);
	mint np = *Args[1].integer;

	vector<numvector<double, 3>> CC;
	CC.resize(np);
	for (size_t i = 0; i < np; ++i)
		CC[i] = { pts[2 * i], pts[2 * i + 1], 0.0 }; //Вектор - трехмерный!!!

	vector<vector<vector<double>>> ME;

	tLayerAver_Matrskos(CC, ME, "mtr");

	mint dimsout[] = { 2 * np+1, 2 * np+1 };

	MTensor v;
	libData->MTensor_new(MType_Real, 2, dimsout, &v);
	
	for (mint p = 0; p < 2 * np + 1; ++p)
	for (mint s = 0; s < 2 * np + 1; ++s)
		((double*)(v->data))[p*(2 * np + 1) + s] = ME[p][s][0];
		
	MArgument_setMTensor(Res, v);
	return LIBRARY_NO_ERROR;
}


//Возвращает матрицу I00
EXTERN_C DLLEXPORT int I00(WolframLibraryData libData,
	mint Argc, MArgument *Args, MArgument Res)
{
	double* pts = (double*)(MArgument_getMTensor(Args[0])->data);
	mint np = *Args[1].integer;

	vector<numvector<double, 3>> CC;
	CC.resize(np);
	for (size_t i = 0; i < np; ++i)
		CC[i] = { pts[2 * i], pts[2 * i + 1], 0.0 }; //Вектор - трехмерный!!!

	vector<vector<vector<double>>> ME;

	tLayerAver_Matrskos(CC, ME, "00");

	mint dimsout[] = { np, np, 3 };

	MTensor v;
	libData->MTensor_new(MType_Real, 3, dimsout, &v);
	
	for (mint p = 0; p < np; ++p)
	for (mint s = 0; s < np; ++s)
	for (mint q = 0; q < 3; ++q)
		((double*)(v->data))[p*(3*np) + s*3 + q] = ME[p][s][q];
		
	MArgument_setMTensor(Res, v);
	return LIBRARY_NO_ERROR;
}

//Возвращает матрицу I01
EXTERN_C DLLEXPORT int I01(WolframLibraryData libData,
	mint Argc, MArgument *Args, MArgument Res)
{
	double* pts = (double*)(MArgument_getMTensor(Args[0])->data);
	mint np = *Args[1].integer;

	vector<numvector<double, 3>> CC;
	CC.resize(np);
	for (size_t i = 0; i < np; ++i)
		CC[i] = { pts[2 * i], pts[2 * i + 1], 0.0 }; //Вектор - трехмерный!!!

	vector<vector<vector<double>>> ME;

	tLayerAver_Matrskos(CC, ME, "01");

	mint dimsout[] = { np, np, 3 };

	MTensor v;
	libData->MTensor_new(MType_Real, 3, dimsout, &v);
	
	for (mint p = 0; p < np; ++p)
	for (mint s = 0; s < np; ++s)
	for (mint q = 0; q < 3; ++q)
		((double*)(v->data))[p*(3*np) + s*3 + q] = ME[p][s][q];
		
	MArgument_setMTensor(Res, v);
	return LIBRARY_NO_ERROR;
}

//Возвращает матрицу I10
EXTERN_C DLLEXPORT int I10(WolframLibraryData libData,
	mint Argc, MArgument *Args, MArgument Res)
{
	double* pts = (double*)(MArgument_getMTensor(Args[0])->data);
	mint np = *Args[1].integer;

	vector<numvector<double, 3>> CC;
	CC.resize(np);
	for (size_t i = 0; i < np; ++i)
		CC[i] = { pts[2 * i], pts[2 * i + 1], 0.0 }; //Вектор - трехмерный!!!

	vector<vector<vector<double>>> ME;

	tLayerAver_Matrskos(CC, ME, "10");

	mint dimsout[] = { np, np, 3 };

	MTensor v;
	libData->MTensor_new(MType_Real, 3, dimsout, &v);
	
	for (mint p = 0; p < np; ++p)
	for (mint s = 0; s < np; ++s)
	for (mint q = 0; q < 3; ++q)
		((double*)(v->data))[p*(3*np) + s*3 + q] = ME[p][s][q];
		
	MArgument_setMTensor(Res, v);
	return LIBRARY_NO_ERROR;
}

//Возвращает матрицу I11
EXTERN_C DLLEXPORT int I11(WolframLibraryData libData,
	mint Argc, MArgument *Args, MArgument Res)
{
	double* pts = (double*)(MArgument_getMTensor(Args[0])->data);
	mint np = *Args[1].integer;

	vector<numvector<double, 3>> CC;
	CC.resize(np);
	for (size_t i = 0; i < np; ++i)
		CC[i] = { pts[2 * i], pts[2 * i + 1], 0.0 }; //Вектор - трехмерный!!!

	vector<vector<vector<double>>> ME;

	tLayerAver_Matrskos(CC, ME, "11");

	mint dimsout[] = { np, np, 3 };

	MTensor v;
	libData->MTensor_new(MType_Real, 3, dimsout, &v);
	
	for (mint p = 0; p < np; ++p)
	for (mint s = 0; s < np; ++s)
	for (mint q = 0; q < 3; ++q)
		((double*)(v->data))[p*(3*np) + s*3 + q] = ME[p][s][q];
		
	MArgument_setMTensor(Res, v);
	return LIBRARY_NO_ERROR;
}