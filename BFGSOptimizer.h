/*
@ author: luchang li, luchangli1993@163.com
@ huazhong university of science and technology
@ 2017/12/07
@ free use for research
*/

#pragma once

#include <math.h>

#define Max(a, b)    (((a) > (b)) ? (a) : (b))
#define Min(a, b)    (((a) < (b)) ? (a) : (b))

/*
IType: type of input data for optimize
ParaNum: number of optimize parameter

IterateNum: total iteration number
  2 ParaNum: 4

  3 ParaNum: 6
  4-5 ParaNum: 8
  6-7 ParaNum: 11

IterateNum_bs: bisection iteration to find best walk length
  IterateNum_bs = 11
*/

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
class BFGSOptimizer {
 public:

  float FitPara[ParaNum]; // cur value

 private :
  float grad[ParaNum];  // gradient
  float d0[ParaNum];  // direction
  float D0[ParaNum * ParaNum]; // inv of Hessian matrix

  float sk[ParaNum]; // bfgs quasi-newton method
  float yk[ParaNum];

  float Para_L[ParaNum];
  float Para_U[ParaNum];
  int ConstrainType; // 0x01|0x02

  // function pointer to get users'  function  for specific curve fitting
  void (* pPreFitting)(float* FitPara, IType* ix, IType* iy, int DataNum);
  float (* pTargerF)(float* FitPara, IType* ix, IType* iy, int DataNum);

  float ScalingCoeff[ParaNum];

 public:

  // function pointer to get users'  function  for specific curve fitting
  // if ipPreFitting is NULL, then you must pre set FitPara manually after the class is created
  // ipTargerF must not be NULL
  BFGSOptimizer(void(* ipPreFitting)(float* FitPara, IType* ix, IType* iy, int DataNum),
                float(* ipTargerF)(float* FitPara,
                                   IType* ix, IType* iy, int DataNum)) {
    int cnt;

    pPreFitting = ipPreFitting;
    pTargerF = ipTargerF;

    ConstrainType = 0;
    for (cnt = 0; cnt < ParaNum; cnt++) {
      ScalingCoeff[cnt] = 1.0f; // without scaling
      FitPara[0] = 0.0f;
    }
  }

  void BFGSOptimize(IType* ix, IType* iy, int DataNum);

  // set parameters value range constrains
  void SetParaConstrain(float* iPara_L, float iPara_U, int iConstrainType);

  void PrintfFitPara(char* pstr) {
    int cnt;
    printf("%s fit result:", pstr);
    for (cnt = 0; cnt < ParaNum; cnt++) {
      printf("%f ", FitPara[cnt]);
    }
    printf("\n");
  }

 private:

  void BFGSOptimizer_Core(IType* ix, IType* iy, int DataNum);
  void PreFitting(float* FitPara, IType* ix, IType* iy, int DataNum);
  float TargerFunction(float* FitPara, IType* ix, IType* iy, int DataNum);
  void GradientCalc(float* FitPara, float* grad, IType* ix, IType* iy, int DataNum);

  // scaling can enlarge parameters' range
  void GetScalingCoeff(float* iFitPara, float* ScalingCoeff);
  void ApplyScalingCoeff(float* iFitPara, float* ScalingCoeff);
  void ApplyRScalingCoeff(float* iFitPara, float* ScalingCoeff);

  void D0Init(float* D0);
  void UpdatePara(float* oX0, float* iX0, float* d0, float coeff);
  void ConstrainPara();
  void UpdateInvHessian(float* D0, float* sk, float* yk);
  void MatMultiplyVector(float* D0, float* grad, float* d0);
};

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::BFGSOptimize(IType* ix, IType* iy, int DataNum) {
  int rcnt;
  // initial guess

  if (pPreFitting != NULL) {
    PreFitting(FitPara, ix, iy, DataNum);

  } else {
//    printf("null pre fit");
  }

//  PrintfFitPara("initial");

  GetScalingCoeff(FitPara, ScalingCoeff);
  ApplyScalingCoeff(FitPara, ScalingCoeff);

  // gradient calculation
  GradientCalc(FitPara, grad, ix, iy, DataNum);

  // first search direction
  for (rcnt = 0; rcnt < ParaNum; rcnt++) {
    d0[rcnt] = -grad[rcnt];
  }

  D0Init(D0);

  BFGSOptimizer_Core(ix, iy, DataNum);
  ApplyRScalingCoeff(FitPara, ScalingCoeff);
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::BFGSOptimizer_Core(IType* ix, IType* iy, int DataNum) {
  // adjust d0
  float td0_total;

  float tgrad[ParaNum];

  int itcnt = 0; // iteration number
  int bcnt = 0;
  int rcnt;

  float scale;
  float xd[ParaNum * 2];

  float ddat[2];
  float dpos[2];

  int xdsel = 0;

  for (itcnt = 0; itcnt < IterateNum; itcnt++) {
    // adjust d0
    td0_total = 0;
    for (rcnt = 0; rcnt < ParaNum; rcnt++) {
      td0_total += abs(d0[rcnt]);
    }

    td0_total = td0_total / ParaNum;

    // normalize
    for (rcnt = 0; rcnt < ParaNum; rcnt++) {
      d0[rcnt] = d0[rcnt] / td0_total;
    }

    dpos[0] = 0.00001f; // scale factor left limit, should not equal to 0 and smaller
    dpos[1] = 1.0f; // scale factor right limit, should not lager than 2

    UpdatePara(&xd[0], FitPara, d0, 0.0001f);
    UpdatePara(&xd[ParaNum], FitPara, d0, 1.0f);

    ddat[0] = TargerFunction(&xd[0], ix, iy, DataNum);
    ddat[1] = TargerFunction(&xd[ParaNum], ix, iy, DataNum);

    // bisection method to find best walk length
    for (bcnt = 0; bcnt < IterateNum_bs; bcnt++) {
      // which part shrink
      xdsel = (ddat[0] < ddat[1]);

      dpos[xdsel] = (dpos[0] + dpos[1]) * 0.5f; //  /2.0f which one shift

      if (bcnt < IterateNum_bs - 1) {
        UpdatePara(&xd[xdsel * ParaNum], FitPara, d0, dpos[xdsel]); // xd=ininf+d0*scale
        ddat[xdsel] = TargerFunction(&xd[xdsel * ParaNum], ix, iy, DataNum);
      }
    }

    scale = (dpos[0] + dpos[1]) * 0.5f; //  /2.0f calculated direction


    for (rcnt = 0; rcnt < ParaNum; rcnt++) {
      sk[rcnt] = d0[rcnt] * scale;
      FitPara[rcnt] = FitPara[rcnt] + sk[rcnt];
    }

    ConstrainPara(); //

    if (itcnt < IterateNum - 1) {
      for (rcnt = 0; rcnt < ParaNum; rcnt++) {
        tgrad[rcnt] = grad[rcnt];
      }

      GradientCalc(FitPara, grad, ix, iy, DataNum);

      for (rcnt = 0; rcnt < ParaNum; rcnt++) {
        yk[rcnt] = grad[rcnt] - tgrad[rcnt];
      }

      UpdateInvHessian(D0, sk, yk);
      MatMultiplyVector(D0, grad, d0);
    }
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::PreFitting(float* FitPara, IType* ix, IType* iy,
                                                                          int DataNum) {

  pPreFitting(FitPara, ix, iy, DataNum);
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
float BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::TargerFunction(float* FitPara, IType* ix, IType* iy,
                                                                               int DataNum) {
  float tFitPara[ParaNum];
  memcpy(tFitPara, FitPara, ParaNum * sizeof(float));

  ApplyRScalingCoeff(tFitPara, ScalingCoeff);

  return pTargerF(tFitPara, ix, iy, DataNum);
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::GetScalingCoeff(float* iFitPara, float* ScalingCoeff) {
  int cnt = 0;
  for (cnt = 0; cnt < ParaNum; cnt++) {
    ScalingCoeff[cnt] = Max(abs(iFitPara[cnt]), 1.0f);
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::ApplyScalingCoeff(float* iFitPara, float* ScalingCoeff) {
  int cnt = 0;
  for (cnt = 0; cnt < ParaNum; cnt++) {
    iFitPara[cnt] = iFitPara[cnt] / ScalingCoeff[cnt];
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::ApplyRScalingCoeff(float* iFitPara,
                                                                                  float* ScalingCoeff) {
  int cnt = 0;
  for (cnt = 0; cnt < ParaNum; cnt++) {
    iFitPara[cnt] = iFitPara[cnt] * ScalingCoeff[cnt];
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::GradientCalc(float* FitPara, float* grad, IType* ix,
                                                                            IType* iy, int DataNum) {
  float tx0[ParaNum];
  float tgradn;
  float tgradp;
  int cnt;

  for (cnt = 0; cnt < ParaNum; cnt++) {
    tx0[cnt] = FitPara[cnt];
  }

  for (cnt = 0; cnt < ParaNum; cnt++) {
    tx0[cnt] = tx0[cnt] - 0.002f;
    tgradn = TargerFunction(tx0, ix, iy, DataNum);

    tx0[cnt] = tx0[cnt] + 0.004f;
    tgradp = TargerFunction(tx0, ix, iy, DataNum);

    grad[cnt] = (tgradp - tgradn) * 250.0f; // /0.004f;
    tx0[cnt] = tx0[cnt] - 0.002f;
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::D0Init(float* D0) {
  float(* pD0)[ParaNum] = (float (*)[ParaNum]) &D0[0];

  memset(D0, 0, 2 * ParaNum * sizeof(float));

  int rcnt;
  //initial D0
  for (rcnt = 0; rcnt < ParaNum; rcnt++) {
    pD0[rcnt][rcnt] = 1.0f;
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::UpdatePara(float* oX0, float* iX0, float* d0,
                                                                          float coeff) {
  int cnt = 0;
  // oX0=d0*coeff + iX0;
  for (cnt = 0; cnt < ParaNum; cnt++) {
    oX0[cnt] = d0[cnt] * coeff + iX0[cnt];
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::SetParaConstrain(float* iPara_L, float iPara_U,
                                                                                int iConstrainType) {
  ConstrainType = iConstrainType;

  // constrain their lower limit
  if (ConstrainType & 0x01) {
    memcpy(Para_L, iPara_L, ParaNum);
  }

  // constrain their upper limit
  if (ConstrainType & 0x02) {
    memcpy(Para_U, iPara_U, ParaNum);
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::ConstrainPara() {
  int cnt = 0;

  // constrain their lower limit
  if (ConstrainType & 0x01) {
    for (cnt = 0; cnt < ParaNum; cnt++) {
      FitPara[cnt] = Max(FitPara[cnt], Para_L[cnt] / ScalingCoeff[cnt]);
    }
  }

  // constrain their upper limit
  if (ConstrainType & 0x02) {
    for (cnt = 0; cnt < ParaNum; cnt++) {
      FitPara[cnt] = Min(FitPara[cnt], Para_U[cnt] / ScalingCoeff[cnt]);
    }
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::UpdateInvHessian(float* D0, float* sk, float* yk) {
  float divdat;
  float skyk[ParaNum * ParaNum]; //  I-(sk*yk')/(yk'*sk)
  //  float yksk[25]; // the same with skyk but transposition
  float sksk[ParaNum * ParaNum];
  float tD0[ParaNum * ParaNum];

  float(* pskyk)[ParaNum] = (float (*)[ParaNum]) skyk;
  float(* psksk)[ParaNum] = (float (*)[ParaNum]) sksk;
  float(* ptD0)[ParaNum] = (float (*)[ParaNum]) tD0;
  float(* pD0)[ParaNum] = (float (*)[ParaNum]) &D0[0];

  int row = 0;
  int col = 0;
  int cnt = 0;

  // = sk.*yk
  divdat = 0;

  for (cnt = 0; cnt < ParaNum; cnt++) {
    divdat += sk[cnt] * yk[cnt];
  }

  divdat = 1.0f / divdat;
  // divdat = __fdividef(1.0f , divdat); //

  float tdat10[ParaNum];
  float tdat20[ParaNum];

  for (cnt = 0; cnt < ParaNum; cnt++) {
    tdat10[cnt] = yk[cnt] * divdat;
    tdat20[cnt] = sk[cnt] * divdat;
  }

  for (row = 0; row < ParaNum; row++) {

    for (col = 0; col < ParaNum; col++) {
      // I-(sk*yk')/(yk'*sk)
      if (row == col) pskyk[row][col] = 1.0f - sk[row] * tdat10[col];
      else pskyk[row][col] = 0.0f - sk[row] * tdat10[col];

      // (sk*sk')/(yk'*sk)
      psksk[row][col] = sk[row] * tdat20[col];
    }
  }
  //

  for (row = 0; row < ParaNum; row++) {

    for (col = 0; col < ParaNum; col++) {
      // tD0 = skyk*D0
      ptD0[row][col] = 0;

      for (cnt = 0; cnt < ParaNum; cnt++) {
        ptD0[row][col] += pskyk[row][cnt] * pD0[cnt][col];
      }
    }
  }

  for (row = 0; row < ParaNum; row++) {

    for (col = 0; col < ParaNum; col++) {
      // D0 = D0*yksk
      pD0[row][col] = 0;

      for (cnt = 0; cnt < ParaNum; cnt++) {
        pD0[row][col] += ptD0[row][cnt] * pskyk[col][cnt];
      }
    }
  }
  // D0=D0+sksk;

  for (row = 0; row < ParaNum; row++) {

    for (cnt = 0; cnt < ParaNum; cnt++) {
      pD0[row][cnt] = pD0[row][cnt] + psksk[row][cnt];
    }
  }
}

template<class IType, int ParaNum, int IterateNum, int IterateNum_bs>
void BFGSOptimizer<IType, ParaNum, IterateNum, IterateNum_bs>::MatMultiplyVector(float* D0, float* grad, float* d0) {
  // d0=-D0*grad;    %search direction

  int row;
  int cnt = 0;

  float(* pD0)[ParaNum] = (float (*)[ParaNum]) &D0[0];

  for (row = 0; row < ParaNum; row++) {
    d0[row] = 0;

    for (cnt = 0; cnt < ParaNum; cnt++) {
      d0[row] -= pD0[row][cnt] * grad[cnt];
    }
  }
}

