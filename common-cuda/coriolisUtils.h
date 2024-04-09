
#ifndef CORIOLIS_UTILS_H
#define CORIOLIS_UTILS_H

/**
 * References
 * E. M. Goldfield, S. K. Gray  Computer Physics Communications, 98, 1-14 (1996)
 * A. J. H.M. Meijer, E. M. Goldfield  J. Chem. Phys., 108, 5404-5413 (1998)
 * H. Wang, W. H. Thompson, W. H. Miller  J. Phys. Chem. A102, 9372-9379 (1998)
 * H. Wang, W. H. Thompson, W. H. Miller  J. Chem. Phys. 107, 7194 (1998)
 * C. Leforestier J. Chem. Phys. 94, 6388-6397 (1991)
 * John Zhang, Theory and Application of Quantum Molecular Dyanmics, P343
 **/

namespace coriolisUtils {
  
  __device__ __host__ inline int kronecker_delta(const int a, const int b)
  { return a == b ? 1 : 0; }
  
  // sign should be +1 or -1
  __device__ __host__ inline double lambda(const int a, const int b, const int sign)
  { return a >= b ? sqrt(double(a*(a+1)-b*(b+sign))) : 0; }
  
  __device__ __host__ inline double coriolis(const int J, const int l, const int Omega, const int Omega1)
  {
    if(Omega1 == Omega+1) 
      return -sqrt(double(1 + kronecker_delta(Omega, 0)))*lambda(J, Omega, 1)*lambda(l, Omega, 1);
    else if(Omega1 == Omega-1) 
      return -sqrt(double(1 + kronecker_delta(Omega, 1)))*lambda(J, Omega, -1)*lambda(l, Omega, -1);
    else
      return 0;
  }
};

#endif /* CORIOLIS_UTILS_H */
