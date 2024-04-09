
#ifndef CUDA_MATH_H
#define CUDA_MATH_H

#if defined __NVCC__ || defined __HIPCC__

#include <cuda_runtime.h>
#include <cufft.h>

#include "complex.h"

namespace cudaMath {

  inline bool is_pow_2(int x) { return ((x&(x-1)) == 0); }
  
  __device__ __host__ inline double sq(const double x) { return x*x; }

#if 0
  __device__ inline double atomicAdd(double *address, double val)
  {
    unsigned long long int* address_as_ull = (unsigned long long int *) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
		      __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
  }
  
  __device__ inline double atomicAdd(double &address, const double &val)
  { return atomicAdd(&address, val); }
  
  __device__ inline Complex atomicAdd(Complex &sum, const Complex &val)
  { 
    atomicAdd(sum.re, val.re); 
    atomicAdd(sum.im, val.im); 
    return sum;
  }
  
  __device__ inline Complex atomicAdd(Complex *sum, const Complex &val)
  { return atomicAdd(*sum, val); }
#endif
  
  template<class T1, class T2, class T3> 
  __global__ void _vector_multiplication_(T1 *vOut, const T2 *vIn1, const T3 *vIn2, const int n)
  {
    const int index = threadIdx.x + blockDim.x*blockIdx.x;
    if(index < n) vOut[index] = vIn1[index]*vIn2[index];
  }
  
  template __global__ 
  void _vector_multiplication_<Complex, Complex, double>(Complex *vOut, 
							 const Complex *vIn1, const double *vIn2, const int n);

  __device__ inline void setup_momentum_for_fft(double *p, const int n, const double xl)
  {
    if(n/2*2 != n) return;
    
    const double two_pi_xl = 2.0*Pi/xl;
    
    for(int i = threadIdx.x; i < n; i += blockDim.x) {
      if(i <= n/2) 
	p[i] = two_pi_xl*i;
      else if(i > n/2)
	p[i] = two_pi_xl*(-n+i);
    }
  }
  
  __device__ inline void setup_kinetic_energy_for_fft(double *kin, const int n, const double xl, 
						      const double mass)
  {
    if(n/2*2 != n) return;
    
    const double two_pi_xl = 2.0*Pi/xl;
    const double mass_inv_half = 0.5/mass;
    
    for(int i = threadIdx.x; i < n; i += blockDim.x) {
      if(i <= n/2) {
	kin[i] = sq(two_pi_xl*i)*mass_inv_half;
      } else if(i > n/2) {
	kin[i] = sq(two_pi_xl*(-n+i))*mass_inv_half;
      }
    }
  }

  __device__ inline void setup_kinetic_energy_for_fft_nonnegative(double *kin, const int n, const double xl, 
								  const double mass)
  {
    if(n/2*2 != n) return;
    
    const double two_pi_xl = 2.0*Pi/xl;
    const double mass_inv_half = 0.5/mass;
    
    for(int i = threadIdx.x; i <= n/2; i += blockDim.x) {
      kin[i] = sq(two_pi_xl*i)*mass_inv_half;
    }
  }
  
  __device__ inline void setup_moments_of_inertia(double *I, const int n, const double r_left, 
						  const double dr, const double mass)
  {
    const double two_mass = mass + mass;
    
    for(int i = threadIdx.x; i < n; i += blockDim.x) {
      const double r = r_left + i*dr;
      I[i] = 1.0/(two_mass*r*r);
    }
  }

  __device__ __host__ inline double WoodsSaxon(const double x, const double Cd, const double xd)
  {
    return 1.0/(1.0 + exp(Cd*(x-xd)));
  }
}

#endif /* __NVCC__ */
#endif /* CUMATH_H */
