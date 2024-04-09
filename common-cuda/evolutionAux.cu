
#ifndef EVOLUTION_AUX_H
#define EVOLUTION_AUX_H

#include "evolutionUtils.h"
#include "cudaMath.h"
#include "coriolisUtils.h"

#ifdef printf
#undef printf
#endif

static __global__ void _print_constant_memory_()
{
#if defined __NVCC__
  printf(" %f %f %f %d\n", r1_dev.left, r1_dev.dr, r1_dev.mass, r1_dev.n);
  printf(" %f %f %f %d\n", r2_dev.left, r2_dev.dr, r2_dev.mass, r2_dev.n);
  for(int i = 0; i < 500; i+=10) 
    printf("%d %18.15f %18.15f\n", i+1, r1_dev.dump[i], r2_dev.dump[i]);
#endif
}

static __global__ void _print_gradient_coeffients_(const int n)
{
#if defined __NVCC__
  for(int i = 0; i < n; i++)
    printf(" %d %18.15e\n", i, gradient_coeffients_dev[i]);
#endif
}

static __global__ void _print_energies_(const int n)
{
#if defined __NVCC__
  for(int i = 0; i < n; i++)
    printf(" %d %18.15f\n", i, energies_dev[i]);
#endif
}

static __global__ void _psi_times_kinetic_energy_(
  Complex *psi_out, const Complex *psi_in, 
	const int n1, const int n2, const int n_theta)
{
  extern __shared__ double kinetic_data[];
  
  double *kin1 = (double *) kinetic_data;
  double *kin2 = &kin1[n1/2+1];
  
  cudaMath::setup_kinetic_energy_for_fft_nonnegative(kin1, r1_dev.n, r1_dev.n*r1_dev.dr, r1_dev.mass);
  cudaMath::setup_kinetic_energy_for_fft(kin2, r2_dev.n, r2_dev.n*r2_dev.dr, r2_dev.mass);

  __syncthreads();
  
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < (n1/2+1)*n2*n_theta) {
    int i = -1; int j = -1; int k = -1;
    cudaUtils::index_2_ijk(index, n1/2+1, n2, n_theta, i, j, k);

#if 0
    const double e = kin1[i] + kin2[j];
    if(e <= kinetic_cutoff) {
      psi_out[index] = e*psi_in[index];
    } else {
      psi_out[index].zero();
    }
#endif
    
    psi_out[index] = (kin1[i] + kin2[j])*psi_in[index];
  }
}

static __global__ void _add_T_radial_weighted_psi_to_H_weighted_psi_(double *HPsi, const double *TRadPsi,
								     const int n1, const int n2, 
								     const int n_theta)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < (n1/2+1)*2*n2*n_theta) {
    int i = -1; int j = -1; int k = -1;
    cudaUtils::index_2_ijk(index, (n1/2+1)*2, n2, n_theta, i, j, k);
    if(i < n1) {
      const int index2 = cudaUtils::ijk_2_index(n1, n2, n_theta, i, j, k);
      HPsi[index2] += TRadPsi[index]/(n1*n2);
    }
  }
}

static __global__ void _add_potential_weighted_psi_to_H_weighted_psi_(double *HPsi, const double *psi,
								      const double *pot, const int n)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n && pot[index] <= potential_cutoff)
    HPsi[index] += pot[index]*psi[index];
}

static __global__ void _add_T_bend_T_sym_to_T_angle_legendre_psi_dev_(double *TangPsi, const double *psi,
								      const int n1, const int n2, 
								      const int nLegs,
								      const int J, const int omega,
								      const int a, const int b)
{
  extern __shared__ double rotational_moments[];

  double *I1 = rotational_moments;
  double *I2 = &I1[n1];
  double &Tsym = I2[n2];
  
  cudaMath::setup_moments_of_inertia(I1, r1_dev.n, r1_dev.left, r1_dev.dr, r1_dev.mass);
  cudaMath::setup_moments_of_inertia(I2, r2_dev.n, r2_dev.left, r2_dev.dr, r2_dev.mass);

  if(threadIdx.x == 0) Tsym = double(J*(J+1) - 2*omega*omega);

  __syncthreads();
  
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n1*n2*nLegs) {
    int i = -1; int j = -1; int l = -1;
    cudaUtils::index_2_ijk(index, n1, n2, nLegs, i, j, l);
    // l += omega;
    l = a*l + b;
    
    const double e = (I1[i]+I2[j])*l*(l+1) + I1[i]*Tsym;
    if(e <= potential_cutoff) TangPsi[index] += e*psi[index];
    //TangPsi[index] += ((I1[i]+I2[j])*l*(l+1) + I1[i]*Tsym)*psi[index];
  }
}

static __global__ void _add_T_asym_to_T_angle_legendre_psi_dev_(double *TangPsi, const double *psi,
								const int n1, const int n2, 
								const int nLegs,
								const int J,
								const int Omega, const int Omega1,
								const int OmegaMax,
								const int a, const int b)
{
  extern __shared__ double I1[];
  
  cudaMath::setup_moments_of_inertia(I1, r1_dev.n, r1_dev.left, r1_dev.dr, r1_dev.mass);
  
  __syncthreads();
  
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n1*n2*nLegs) {
    int i = -1; int j = -1; int l = -1;
    cudaUtils::index_2_ijk(index, n1, n2, nLegs, i, j, l);
    //l += OmegaMax;
    l = a*l + b;
    const double c = coriolisUtils::coriolis(J, l, Omega, Omega1);
    const double e = I1[i]*c;
    if(e <= potential_cutoff) TangPsi[index] += e*psi[index];
    //TangPsi[index] += I1[i]*c*psi[index];
  }
}

static __global__ void _dump_wavepacket_(double *psi, const int n1, const int n2, const int n_theta)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  
  if(index < n1*n2*n_theta) {
    int i = -1; int j = -1; int k = -1;
    cudaUtils::index_2_ijk(index, n1, n2, n_theta, i, j, k);
    psi[index] *= r1_dev.dump[i]*r2_dev.dump[j];
  }
}

static __global__ void _daxpy_(double *y, const double *x, const double alpha, const double beta,
			       const int n)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n)
    y[index] = alpha*x[index] + beta*y[index];
}

#if 0
static __global__ void _setup_potential_scale_(int *scale, const double *pot_dev,
					       const double cutoff, const int n)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n) 
    scale[index] = pot_dev[index] < cutoff ? 1 : 0;
}

static __global__ void _scale_wavepacket_with_potential_cutoff_(double *psi, const double *potential,
								const double cutoff, const int n)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n && potential[index] > cutoff) psi[index] = 0.0;
}
#endif

static __global__ void _psi_time_to_fai_energy_on_dividing_surface_
(const int n, const int n_energies,
 const double t, const double dt,
 const double *psi_real_dev, const double *psi_imag_dev, 
 const double *d_psi_real_dev, const double *d_psi_imag_dev,
 Complex *fai_dev, Complex *d_fai_dev)
{
  extern __shared__ Complex expIEtDt[];
  
  for(int i = threadIdx.x; i < n_energies; i += blockDim.x)
    expIEtDt[i] = exp(Complex(0.0, t*energies_dev[i]))*dt;
  
  __syncthreads();

  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  if(index < n*n_energies) {
    int i = -1; int iE = -1;
    cudaUtils::index_2_ij(index, n, n_energies, i, iE);
    fai_dev[index] += expIEtDt[iE]*Complex(psi_real_dev[i], psi_imag_dev[i]);
    d_fai_dev[index] += expIEtDt[iE]*Complex(d_psi_real_dev[i], d_psi_imag_dev[i]);
  }
}

#endif /* EVOLUTION_AUX_H */

