
/* $Id$ */

#include <iostream>
#include <cassert>
#include <cstring>

#include <mex.h>
#include "fort.h"
#include "matlabUtils.h"

extern "C" void FORT(bkmp2)(const double *r, double &v, double *dv, const int &id);


inline void bkmp2(const double &r1, const double &r2, const double &r3,
		  double &v, double *dv, const int &id)
{
  double r[] = { r1, r2, r3 };
  FORT(bkmp2)(r, v, dv, id);
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 3) 
    mexErrMsgTxt("BKMP2Mex requires 3 input arguments");
  
  if (nlhs > 3)
    mexErrMsgTxt("BKMP2Mex output arguments <= 2");
  
  int i = 0;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  const int m = mxGetM(prhs[i]);
  const int n = mxGetN(prhs[i]);
  const double *r1 = mxGetPr(prhs[i]);
  insist(r1);

  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *r2 = mxGetPr(prhs[i]);
  insist(r2);
  
  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *r3 = mxGetPr(prhs[i]);
  insist(r3);

  if(nlhs == 0 || nlhs == 1) {
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *pl0 = mxGetPr(plhs[0]);
    insist(pl0);
    
#pragma omp parallel for if(m*n > 100)		\
  default(shared) schedule(static, 1)		      
    for(int i = 0; i < m*n; i++) {
      double &v = pl0[i];
      double dv[] = { 0.0, 0.0, 0.0 };
      bkmp2(r1[i], r2[i], r3[i], v, dv, 0);
    }
  }

 /*
  else if(nlhs == 2) {
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *pl0 = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(3*m, n, mxREAL);
    double *pl1 = mxGetPr(plhs[1]);
    
#pragma omp parallel for if(m*n > 100)		\
  default(shared) schedule(static, 1)		      
    for(int i = 0; i < m*n; i++) {
      double &v = pl0[i];
      double *dv = pl1+3*i;
      bkmp2(r1[i], r2[i], r3[i], v, dv, 1);
    }
  }
  */
  
  std::cout.flush();
}
