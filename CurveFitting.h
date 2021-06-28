/*
@ author: luchang li, luchangli1993@163.com
@ huazhong university of science and technology
@ 2017/12/07
@ free use for research
*/

#include "BFGSOptimizer.h"

// initial guess function for optimization/ curve fitting,DataNum is the array length
// void xx_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
// optimization target function
// float xx_TargerF(float *FitPara, float *ix, float *iy, int DataNum);
// x can be time,count or just NULL, depend on how you define them

/*
exp fitting
y=a*exp(b*x)
FitPara=[a,b]
*/

void ExpFit_PreFitting(float* FitPara, float* ix, float* iy, int DataNum);
float ExpFit_TargerF(float* FitPara, float* ix, float* iy, int DataNum);

/*
Gaussian Fitting 1 0
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/

void GausFit10_PreFitting(float* FitPara, float* ix, float* iy, int DataNum);
float GausFit10_TargerF(float* FitPara, float* ix, float* iy, int DataNum);

/*
Gaussian Fitting 1 1
Y=a*exp(-(x-x0)^2/(2*sigma^2))+b
FitPara=[A,x0,sigma,b]
*/
void GausFit11_PreFitting(float* FitPara, float* ix, float* iy, int DataNum);
float GausFit11_TargerF(float* FitPara, float* ix, float* iy, int DataNum);
