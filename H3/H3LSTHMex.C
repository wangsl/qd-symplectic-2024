
/* $Id$ */

#include <iostream>
#include <cassert>
#include <cstring>

#include <mex.h>
#include "fort.h"
#include "matlabUtils.h"

extern "C" void FORT(h3slth)(const double *x, double &v);

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 3) 
    mexErrMsgTxt("H3LSTHMex requires 3 input arguments");
  
  if (nlhs > 3)
    mexErrMsgTxt("H3LSTHMex output arguments <= 2");
  
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

  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  double *pl0 = mxGetPr(plhs[0]);
  insist(pl0);
    
#pragma omp parallel for if(m*n > 100)		\
  default(shared) schedule(static, 1)		      
  for(int i = 0; i < m*n; i++) {
    double r[] = { r1[i], r2[i], r3[i] };
    double &v = pl0[i];
    FORT(h3slth)(r, v);
  }

  std::cout.flush();
}
