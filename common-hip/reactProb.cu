
#include <iostream>

#include "cudaUtils.h"
#include "matlabData.h"
#include "evolutionUtils.h"
#include "reactProb.h"

#include "evolutionAux.cu"

#ifdef __HIPCC__
#define hipDoubleComplex hipblasDoubleComplex 
#endif

ReactionProbabilities::ReactionProbabilities(const OmegaWavepacket *wavepacket_) :
  wavepacket(wavepacket_),
  fai_on_surface_dev(0),
  d_fai_on_surface_dev(0),
  psi_real_dev(0), d_psi_real_dev(0),
  psi_imag_dev(0), d_psi_imag_dev(0)
{
  reaction_probabilities.resize(MatlabData::crp_parameters()->n_energies, 0);
  setup_data_on_device();
} 

ReactionProbabilities::~ReactionProbabilities()
{ 
  std::cout << " Destroy reaction probabilities data for Omega: "
	    << wavepacket->omega << std::endl;
  
  reaction_probabilities.resize(0);
  
  wavepacket = 0;
  
  psi_real_dev = 0;
  d_psi_real_dev = 0;
  psi_imag_dev = 0;
  d_psi_imag_dev = 0;
  
  _CUDA_FREE_(fai_on_surface_dev);
  _CUDA_FREE_(d_fai_on_surface_dev);
}

void ReactionProbabilities::setup_data_on_device()
{
  std::cout << " Setup reaction probabilities data on device for Omega: " 
	    << wavepacket->omega << std::endl;
  
  const int &n1 = MatlabData::r1()->n;
  const int &n_theta = MatlabData::theta()->n;
  const int &n_energies = MatlabData::crp_parameters()->n_energies;
  
  if(!fai_on_surface_dev) {
    checkCudaErrors(hipMalloc(&fai_on_surface_dev, n1*n_theta*n_energies*sizeof(Complex)));
    checkCudaErrors(hipMemset(fai_on_surface_dev, 0, n1*n_theta*n_energies*sizeof(Complex)));
  }
  
  if(!d_fai_on_surface_dev) {
    checkCudaErrors(hipMalloc(&d_fai_on_surface_dev, n1*n_theta*n_energies*sizeof(Complex)));
    checkCudaErrors(hipMemset(d_fai_on_surface_dev, 0, n1*n_theta*n_energies*sizeof(Complex)));
  }
}

void ReactionProbabilities::calculate_psi_gradients_on_dividing_surface()
{
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  const double &dr2 = MatlabData::r2()->dr;
  
  const int &n_dividing_surface = MatlabData::crp_parameters()->n_dividing_surface;
  const int &n_gradient_points = MatlabData::crp_parameters()->n_gradient_points;
  
  insist(n_dividing_surface < n2);
  
  psi_real_dev = wavepacket->work_dev;
  d_psi_real_dev = psi_real_dev + n1*n_theta;
  psi_imag_dev = d_psi_real_dev + n1*n_theta;
  d_psi_imag_dev = psi_imag_dev + n1*n_theta;
  
  Num1ststGradient::gradients_2_3d(n1, n2, n_theta, 
				   n_dividing_surface, n_gradient_points, dr2, 
				   wavepacket->weighted_psi_real_dev, psi_real_dev, d_psi_real_dev);

  Num1ststGradient::gradients_2_3d(n1, n2, n_theta, 
				   n_dividing_surface, n_gradient_points, dr2, 
				   wavepacket->weighted_psi_imag_dev, psi_imag_dev, d_psi_imag_dev);
}

void ReactionProbabilities::psi_time_to_fai_energy_on_dividing_surface() const
{
  const int &n1 = MatlabData::r1()->n;
  const int &n_theta = MatlabData::theta()->n;
  const int &n_energies = MatlabData::crp_parameters()->n_energies;

  const double &dt = MatlabData::time()->time_step;
  const int &n_steps = MatlabData::time()->steps;
  
  const double t = n_steps*dt;
  
  const int n_threads = _NTHREADS_;
  const int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n_theta*n_energies);
  
  _psi_time_to_fai_energy_on_dividing_surface_<<<n_blocks, n_threads, n_energies*sizeof(Complex)>>>
    (n1*n_theta, n_energies, t, dt, 
     psi_real_dev, psi_imag_dev, d_psi_real_dev, d_psi_imag_dev, 
     fai_on_surface_dev, d_fai_on_surface_dev);
}

void ReactionProbabilities::calculate_reaction_probabilities()
{
  const int &n1 = MatlabData::r1()->n;
  const int &n_theta = MatlabData::theta()->n;
  const double &dr1 = MatlabData::r1()->dr;
  const double &mu2 = MatlabData::r2()->mass;
  
  const int &n_energies = MatlabData::crp_parameters()->n_energies;
  const double *eta_sq = MatlabData::crp_parameters()->eta_sq;
  
  const double dr1_mu2 = dr1/mu2;
  
  insist(reaction_probabilities.size() == n_energies);
  
  for(int iE = 0; iE < n_energies; iE++) {
    const Complex *fai = fai_on_surface_dev + iE*n1*n_theta;
    const Complex *dfai = d_fai_on_surface_dev + iE*n1*n_theta;
    Complex s(0.0, 0.0);
    insist(hipblasZdotc(wavepacket->cublas_handle, n1*n_theta,
		       (const hipDoubleComplex *) dfai, 1,
		       (const hipDoubleComplex *) fai, 1,
		       (hipDoubleComplex *) &s) == HIPBLAS_STATUS_SUCCESS);
    
    reaction_probabilities[iE] = s.imag()/eta_sq[iE]*dr1_mu2;
  }
}

void ReactionProbabilities::calculate_reaction_probabilities(const int calculate)
{
  calculate_psi_gradients_on_dividing_surface();
  
  psi_time_to_fai_energy_on_dividing_surface();
  
  if(calculate) calculate_reaction_probabilities();
}
