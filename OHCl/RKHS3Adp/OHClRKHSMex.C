
#include <iostream>
#include <mex.h>

#include "fort.h"
#include "matlabUtils.h"

#define ATTPACKED(x) __attribute__((packed)) FORT(x)

extern "C" {

  struct {
    int length;
    char data_dir[512];
  } ATTPACKED(rkhsdatadir);
  
  void FORT(prepot)();
  void FORT(pot)(const double &rOH, const double &rHCl,
		 const double &rOCl, double &v);
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray **prhs)
{
  if(nrhs != 3) 
    mexErrMsgTxt("OHClRKHSMex requires 3 input arguments");
  
  if(nlhs > 1)
    mexErrMsgTxt("OHClRKHSMex output arguments == 1");

  static int has_prepot = 0;
  
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
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  double *pl0 = mxGetPr(plhs[0]);
  memset(pl0, 0, n*m*sizeof(mxREAL));

  if(!has_prepot) {
    const char *data_dir = getenv("OHCL_RKHS_DATA_DIR");
    if(!data_dir) 
      mexErrMsgTxt("OHClRKHSMex: no enviorment variable 'OHCL_RKHS_DATA_DIR' defined");
    
    FORT(rkhsdatadir).length = strlen(data_dir);
    memcpy(FORT(rkhsdatadir).data_dir, data_dir, sizeof(char)*strlen(data_dir));
    
    FORT(prepot)();
    has_prepot = 1;
  }

#pragma omp parallel for                        \
  default(shared) schedule(static, 100)		
  for(int i = 0; i < m*n; i++) {
    FORT(pot)(rOH[i], rHCl[i], rOCl[i], pl0[i]);
  }
  std::cout.flush();
}
