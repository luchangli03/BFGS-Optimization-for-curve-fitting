/*
@ author: luchang li, luchangli1993@163.com
@ huazhong university of science and technology
@ 2017/12/07
@ free use for research
*/

#include <iostream>
using namespace std;

#include "BFGSOptimizer.h"

// initial guess function for optimization/ curve fitting
// void xx_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
// optimization target function
// float xx_TargerF(float *FitPara, float *ix, float *iy, int DataNum);

void ZPlaneFit_PreFitting(float *FitPara, float *DataArrayList[3], int DataNum);
float ZPlaneFit_TargerF(float *FitPara, float *DataArrayList[3], int DataNum);

/*
Line order 1 fitting
y=a*x+b

*/

void LineOrder1_PreFitting(float *FitPara, float *DataArrayList[2], int DataNum);

float LineOrder1_TargerF(float *FitPara, float *DataArrayList[2], int DataNum);

