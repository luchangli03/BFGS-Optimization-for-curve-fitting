#include <iostream>
using namespace std;

#include "CurveFitting.h"

int main() {

#define DatLen   15

  // ParaNum, ItNum, ItNum_bs;
  BFGSOptimizer<float, 2, 4, 11> ExpCurveFit(ExpFit_PreFitting, ExpFit_TargerF);

  float ix[DatLen];
  float iy[DatLen] = {0.561983237858735, 0.274092050681491, 0.0936873888873212, 0.0357818099726275, 0.0158309111945142,
                      0.00702655416654909, 0.00437395942094421, 0.00214465107091458, 0.00146739283799419,
                      0.000649039139882044, 0.000705477325958744, 0.000649039139882044, 0.000225752744306798,
                      0.000310410023421847, 0.000169314558230098};

  int cnt;

  for (cnt = 0; cnt < DatLen; cnt++) {
    ix[cnt] = cnt;
  }

  ExpCurveFit.BFGSOptimize(ix, iy, DatLen);
  ExpCurveFit.PrintfFitPara("final");

  int ix1[DatLen];
  int iy1[DatLen] = {1, 1, 2, 4, 6, 7, 9, 9, 9, 9, 7, 5, 4, 2, 1};

  for (cnt = 0; cnt < DatLen; cnt++) {
    ix1[cnt] = cnt;
  }
  BFGSOptimizer<int, 3, 6, 11> GausCurveFit(NULL, GausFit10i_TargerF);
  GausCurveFit.FitPara[0] = 9;
  GausCurveFit.FitPara[1] = 7;
  GausCurveFit.FitPara[2] = 3;

  GausCurveFit.BFGSOptimize(ix1, iy1, DatLen);

  GausCurveFit.PrintfFitPara("final");

  return 0;
}
