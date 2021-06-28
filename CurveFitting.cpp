/*
@ author: luchang li, luchangli1993@163.com
@ huazhong university of science and technology
@ 2017/12/07
@ free use for research
*/

#include "CurveFitting.h"

/*
exp fitting
y=a*exp(b*x)
FitPara=[a,b]

*/

void ExpFit_PreFitting(float* FitPara, float* ix, float* iy, int DataNum) {
  FitPara[0] = iy[0];
  FitPara[1] = logf(iy[1] / iy[0]);
}

float ExpFit_TargerF(float* FitPara, float* ix, float* iy, int DataNum) {

  float a = FitPara[0];
  float b = FitPara[1];

  float y0;

  int cnt = 0;
  float SquareError = 0;

  for (cnt = 0; cnt < DataNum; cnt++) {
    y0 = a * expf(b * ix[cnt]);

    SquareError += powf(y0 - iy[cnt], 2);
  }
  return SquareError;
}

/*
Gaussian Fitting 1 0
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/


void GausFit10_PreFitting(float* FitPara, float* ix, float* iy, int DataNum) {
  float mdat = iy[0];
  float mpos = 0;

  int cnt = 0;
  for (cnt = 0; cnt < DataNum; cnt++) {
    if (mdat < iy[cnt]) {
      mdat = iy[cnt];
      mpos = cnt;
    }
  }
  FitPara[0] = mdat;
  FitPara[1] = mpos;
  FitPara[2] = mpos / 2.5f;

}

float GausFit10_TargerF(float* FitPara, float* ix, float* iy, int DataNum) {
  float a = FitPara[0];
  float x0 = FitPara[1];
  float sigma = FitPara[2];
  float y0;

  int cnt = 0;
  float SquareError = 0;

  for (cnt = 0; cnt < DataNum; cnt++) {
    y0 = a * expf(-(ix[cnt] - x0) * (ix[cnt] - x0) / (2 * sigma * sigma));

    SquareError += powf(y0 - iy[cnt], 2);
  }

  return SquareError;

}

/*
Gaussian Fitting 1 1
Y=a*exp(-(x-x0)^2/(2*sigma^2))+b
FitPara=[A,x0,sigma,b]
*/


void GausFit11_PreFitting(float* FitPara, float* ix, float* iy, int DataNum) {
  float mdat = iy[0];
  float mpos = 0;

  int cnt = 0;

  for (cnt = 0; cnt < DataNum; cnt++) {
    if (mdat < iy[cnt]) {
      mdat = iy[cnt];
      mpos = cnt;
    }
  }

  float b = Min(iy[0], iy[DataNum - 1]);

  FitPara[0] = mdat - b;
  FitPara[1] = mpos;
  FitPara[2] = mpos / 2.5f;
  FitPara[3] = b;

}

float GausFit11_TargerF(float* FitPara, float* ix, float* iy, int DataNum) {
  float a = FitPara[0];
  float x0 = FitPara[1];
  float sigma = FitPara[2];
  float b = FitPara[3];

  float y0;

  int cnt = 0;
  float SquareError = 0;

  for (cnt = 0; cnt < DataNum; cnt++) {
    y0 = a * expf(-(ix[cnt] - x0) * (ix[cnt] - x0) / (2 * sigma * sigma)) + b;

    SquareError += powf(y0 - iy[cnt], 2);
  }

  return SquareError;

}


