#include <iostream>
using namespace std;

#include "CurveFitting.h"


int main()
{
#define DataLen	10
float x[DataLen] = { 0,1,2,3,4,5,6,7,8,9 };
float y[DataLen] = { 0.601981941401637,1.26297128454014,2.65407909847678,3.68921450314001,4.74815159282371,5.45054159850250,6.08382137799693,7.22897696871682,8.91333736150167,9.15237801896922 };

#define DataArrayNum		2

float *DataArray[DataArrayNum] = { x,y };


BFGSOptimizer_TypeDef< float, DataArrayNum, 2, 5, 11>BFGSOptimizer(LineOrder1_PreFitting, LineOrder1_TargerF);

BFGSOptimizer.BFGSOptimize(DataArray, DataLen);

BFGSOptimizer.PrintfFitPara();

return 0;
}



