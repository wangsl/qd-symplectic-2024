
#include <iostream>
#include <cassert>
#include <cstring>
#include <mex.h>

#include "fort.h"
#include "matlabUtils.h"

extern "C" void FORT(ho2dmbeiv)(const double *x, double &v);

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 3) 
    mexErrMsgTxt("DMBEIVMex requires 3 input arguments");
  
  if(nlhs > 1)
    mexErrMsgTxt("DMBEIVMex output arguments == 1");
  
  int i = 0;
  
  i = 0;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  const int m = mxGetM(prhs[i]);
  const int n = mxGetN(prhs[i]);
  const double *rOO = (const double *) mxGetData(prhs[i]);
  insist(rOO);
  
  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *rOH1 = (const double *) mxGetData(prhs[i]);
  insist(rOH1);
  
  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *rOH2 = (const double *) mxGetData(prhs[i]);
  insist(rOH2);
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  double *pl0 = (double *) mxGetData(plhs[0]);
  
#pragma omp parallel for default(shared) schedule(static, 1)
  for(int i = 0; i < m*n; i++) {
    double &v = pl0[i];
    const double R [] = { rOO[i], rOH1[i], rOH2[i] };
    FORT(ho2dmbeiv)(R, v);
  }
  
  std::cout.flush();
}
