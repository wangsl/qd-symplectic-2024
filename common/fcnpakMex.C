
#include <iostream>
#include <cmath>
#include <cassert>

#include <mex.h>

/***************
   Mathematica code
   $MaxExtraPrecision = 1024;
   P[l_, m_, x_] := Sqrt[(2 l + 1)/2 (l - m)!/(l + m)!] LegendreP[l, m, x] (-1)^m
   m = 5; lmax = 200; \[Theta] = 2/5 Pi;
   N[P[Table[l, {l, m, lmax}], m, Cos[\[Theta]]], 20] // TableForm
   
   P_l_m(-x) = (-1)^(l+m)*P_l_m(x)
****************/

#define FORT(x) x##_

extern "C" void FORT(xdlegf)(const double &DNU1, const int &NUDIFF,
			     const int &MU1, const int &MU2,
			     const double &THETA, const int &ID,
			     double *PQA, int *IPQA);

inline void convert_to_decimal(const int &n, double *p, const int *ip)
{
  for(int i = 0; i < n; i++)
    if(ip[i]) p[i] *= pow(10.0, ip[i]);
}

static void fcnpak_4a(const int &nu, const int &mu1, const int &mu2,
		      const double &theta, double *p)
{
  int *ip = new int [mu2-mu1+1];
  assert(ip);

  const double nu_ = nu;
  
  FORT(xdlegf)(nu_, 0, mu1, mu2, theta, 4, p, ip);

  convert_to_decimal(mu2-mu1+1, p, ip);
  
  if(ip) { delete [] ip; ip = 0; }
}

static void fcnpak_4b(const int &nu1, const int &nu2, const int &mu,
		      const double &theta, double *p)
{
  int *ip = new int [nu2-nu1+1];
  assert(ip);

  const double nu1_ = nu1;
  
  FORT(xdlegf)(nu1_, nu2-nu1, mu, mu, theta, 4, p, ip);

  convert_to_decimal(nu2-nu1+1, p, ip);
  
  if(ip) { delete [] ip; ip = 0; }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray **prhs)
{
  if(nrhs != 5) 
    mexErrMsgTxt("fcnpakMex requires 5 input arguments");
  
  if(nlhs > 1)
    mexErrMsgTxt("fcnpakMex output arguments == 1");

  const double half_pi = 0.5*M_PI;
  
  const int &index = *((int *) mxGetData(prhs[0]));
  
  if(index != 0 && index != 1)
    mexErrMsgTxt("fcnpakMex index should be 0 (4a) or 1 (4b) only");

  if(index == 0) { // fcnpak 4a
    const int &nu = *((int *) mxGetData(prhs[1]));
    const int &mu1 = *((int *) mxGetData(prhs[2]));
    const int &mu2 = *((int *) mxGetData(prhs[3]));

    const int m = mxGetM(prhs[4]);
    const int n = mxGetN(prhs[4]);
    const double *x = (const double *) mxGetData(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(mu2-mu1+1, m*n, mxREAL);

    double *p = mxGetPr(plhs[0]);

#pragma omp parallel for default(shared) schedule(static, 1)
    for(int i = 0; i < m*n; i++) {
      double *p_ = p + i*(mu2-mu1+1);
      fcnpak_4a(nu, mu1, mu2, acos(fabs(x[i])), p_);
      if(x[i] < 0) {
	for(int mu = mu1; mu <= mu2; mu++)
	  if((nu+mu)%2) p_[mu-mu1] = -p_[mu-mu1];
      }
    }
    
  } else if (index == 1) { // fcnpak 4b
    const int &nu1 = *((int *) mxGetData(prhs[1]));
    const int &nu2 = *((int *) mxGetData(prhs[2]));
    const int &mu = *((int *) mxGetData(prhs[3]));
    
    const int m = mxGetM(prhs[4]);
    const int n = mxGetN(prhs[4]);
    const double *x = (const double *) mxGetData(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(nu2-nu1+1, m*n, mxREAL);

    double *p = mxGetPr(plhs[0]);

#pragma omp parallel for default(shared) schedule(static, 1)    
    for(int i = 0; i < m*n; i++) {
      double *p_ = p + i*(nu2-nu1+1);
      fcnpak_4b(nu1, nu2, mu, acos(fabs(x[i])), p_);
      if(x[i] < 0) {
	for(int nu = nu1; nu <= nu2; nu++)
	  if((nu+mu)%2) p_[nu-nu1] = -p_[nu-nu1];
      }
    }
  }

  std::cout.flush();
}

