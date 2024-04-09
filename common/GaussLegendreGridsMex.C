
#include <cmath>
#include <cassert>
#include <iostream>

#include <mex.h>

// Gauss-Legendre nodes and weights with kind = 1

#define FORT(x) x##_

extern "C" void FORT(gaussq)(const int &kind, const int &n,
			     const double &alpha, const double &beta,
			     const int &kpts, const double *endpts,
			     double *b, double *t, double *w);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 1) 
    mexErrMsgTxt("GaussLegederMex requires 1 input argument");
  
  const int &n = *((int *) mxGetData(prhs[0]));

  if(n > 599) std::cout << " n is too big, be care of the accuracy: n = "
			<< n << std::endl;
    
  const int kind = 1;
  const double alpha = 0.0;
  const double beta = 0.0;
  const int kpts = 0;
  double *endpts = 0;
  
  double *b = new double [n];

  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
  
  double *t = mxGetPr(plhs[0]);
  double *w = mxGetPr(plhs[1]);

  FORT(gaussq)(kind, n, alpha, beta, kpts, endpts, b, t, w);

  if(b) { delete [] b; b = 0; }
  
  std::cout.flush();
}

