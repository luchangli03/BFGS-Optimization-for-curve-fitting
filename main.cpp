#include <iostream>
using namespace std;

#include "CurveFitting.h"


int main()
{

#define DatLen   15

	// ParaNum, ItNum, ItNum_bs; 
	BFGSOptimizer< 2, 4, 11> ExpCurveFit(ExpFit_PreFitting, ExpFit_TargerF);

	float ix[DatLen];
	float iy[DatLen] = { 0.561983237858735, 0.274092050681491, 0.0936873888873212, 0.0357818099726275, 0.0158309111945142, 0.00702655416654909, 0.00437395942094421, 0.00214465107091458, 0.00146739283799419, 0.000649039139882044, 0.000705477325958744, 0.000649039139882044, 0.000225752744306798, 0.000310410023421847, 0.000169314558230098 };

	int cnt;


	for (cnt = 0; cnt < DatLen; cnt++)
	{
		ix[cnt] = cnt;
	}

	ExpCurveFit.BFGSOptimizer_Top(ix, iy, DatLen);
	ExpCurveFit.PrintfFitPara("final");



	BFGSOptimizer< 4, 8, 11> GausCurveFit(GausFit11_PreFitting, GausFit11_TargerF);

	float ix1[DatLen];
	float iy1[DatLen] = { 1.03787663582584, 1.83038049957306, 2.97498249283382, 4.45630428499173, 6.15194791898853, 7.82704538241868, 9.17762784414991, 9.91768757606545, 9.87728971930251, 9.06593396336589, 7.66892868588658, 5.97866523316379, 4.29557358210739, 2.84436599826735, 1.73579000492169 };

	for (cnt = 0; cnt < DatLen; cnt++)
	{
		ix1[cnt] = cnt;
	}

	GausCurveFit.BFGSOptimizer_Top(ix1, iy1, DatLen);

	GausCurveFit.PrintfFitPara("final");

	return 0;
}
