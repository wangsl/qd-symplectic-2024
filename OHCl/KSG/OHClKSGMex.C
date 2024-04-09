
#include <iostream>
#include <mex.h>

#include "fort.h"
#include "matlabUtils.h"

extern "C" {
  void FORT(ohclksgpes)(const double &rOH, const double &rOCl, 
			const double &rHCl, double &V);
  // InitKSHParameters
  void FORT(initkshparameters)();
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 3) 
    mexErrMsgTxt("OHClKSGMex requires 3 input arguments");
  
  if(nlhs > 1)
    mexErrMsgTxt("OHClKSGMex output arguments == 1");
  
  int i = 0;

  i = 0;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  const int m = mxGetM(prhs[i]);
  const int n = mxGetN(prhs[i]);
  const double *rOH = mxGetPr(prhs[i]);
  insist(rOH);

  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *rOCl = mxGetPr(prhs[i]);
  insist(rOCl);
  
  i++;
  insist(mxIsDouble(prhs[i]));
  insist(mxGetNumberOfDimensions(prhs[i]) == 2);
  insist(m == mxGetM(prhs[i]) && n == mxGetN(prhs[i]));
  const double *rHCl = mxGetPr(prhs[i]);
  insist(rHCl);
  
  // FORT(ohclksgpes) can be run in OpenMP mode only after 
  // declare KSGPESVar as threadprivate
  if(nlhs == 0 || nlhs == 1) {
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *pl0 = mxGetPr(plhs[0]);
    
    FORT(initkshparameters)();
    
#pragma omp parallel for			\
  default(shared) schedule(static, 1)
    for(int i = 0; i < m*n; i++) {
      double &V = pl0[i];
      FORT(ohclksgpes)(rOH[i], rOCl[i], rHCl[i], V);
    }
  }
  std::cout.flush();
}
