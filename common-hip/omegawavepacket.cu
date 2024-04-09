
#include <iostream>

#include "omegawavepacket.h"
#include "cudaUtils.h"
#include "matlabUtils.h"
#include "matlabData.h"
#include "symplecticUtils.h"

#include "evolutionAux.cu"
#include "reactProb.h"

inline int OmegaWavepacket::minimum_Legendres_order() const
{
  int l_min = -1;
  const int &rot_states = MatlabData::options()->rotational_states;
  if(rot_states == _RotStatesAll_) 
    l_min = omega;
  else if(rot_states == _RotStatesOdd_)
    l_min = EvolutionUtils::int_to_odd_right(omega);
  else if(rot_states == _RotStatesEven_)
    l_min = EvolutionUtils::int_to_even_right(omega);

  return l_min;
}

inline int OmegaWavepacket::number_of_associated_Legendres() const
{
  const int &l_max = MatlabData::wavepacket_parameters()->l_max;
  const int &rot_states = MatlabData::options()->rotational_states;
  
  int n_ass_Legs = -1;
  if(rot_states == _RotStatesAll_) 
    n_ass_Legs = l_max - omega + 1;
  else if(rot_states == _RotStatesOdd_)
    n_ass_Legs = (EvolutionUtils::int_to_odd_left(l_max) - 
		  EvolutionUtils::int_to_odd_right(omega))/2 + 1;
  else if(rot_states == _RotStatesEven_)
    n_ass_Legs = (EvolutionUtils::int_to_even_left(l_max) -
		  EvolutionUtils::int_to_even_right(omega))/2 + 1;

  return n_ass_Legs;
}

OmegaWavepacket::OmegaWavepacket(int omega_,
				 const double *potential_dev_, 
				 hipblasHandle_t &cublas_handle_,
				 hipfftHandle &cufft_plan_D2Z_,
				 hipfftHandle &cufft_plan_Z2D_,
				 hipStream_t * &computation_stream_,
				 double * &device_work_dev_
				 ) :
  omega(omega_), 
  potential_dev(potential_dev_),
  cublas_handle(cublas_handle_), 
  cufft_plan_D2Z(cufft_plan_D2Z_),
  cufft_plan_Z2D(cufft_plan_Z2D_),
  computation_stream(computation_stream_),
  device_work_dev(device_work_dev_),
  weighted_psi_real(0), weighted_psi_imag(0), 
  weighted_psi_real_dev(0), weighted_psi_imag_dev(0),
  work_dev(0),
  weighted_associated_legendres_dev(0),
  weighted_psi_dev(0), 
  legendre_psi_dev(0),
  T_angle_legendre_psi_dev(0),
  H_weighted_psi_dev(0),
  crp(0)
{ 
  insist(potential_dev);
  insist(computation_stream);
  copy_weighted_psi_from_host_to_device();
  copy_weighted_associated_legendres_from_host_to_device();
  setup_work_dev();
  // last one to setup reaction probabilities
  setup_reaction_probabilities();
}

OmegaWavepacket::~OmegaWavepacket()
{
  std::cout << " Destroy OmegaWavepacket, Omega: " << omega << std::endl;

  weighted_psi_real = 0;
  weighted_psi_imag = 0;

  potential_dev = 0;

  weighted_psi_dev = 0;
  legendre_psi_dev = 0;
  T_angle_legendre_psi_dev = 0;
  H_weighted_psi_dev = 0;

  _CUDA_FREE_(weighted_psi_real_dev);
  _CUDA_FREE_(weighted_psi_imag_dev);
  _CUDA_FREE_(weighted_associated_legendres_dev);
  _CUDA_FREE_(work_dev);

  if(crp) { delete crp; crp = 0; }
}

void OmegaWavepacket::setup_weighted_psi_dev(const int part)
{
  if(part == _RealPart_)
    weighted_psi_dev = weighted_psi_real_dev;
  else if(part == _ImagPart_)
    weighted_psi_dev = weighted_psi_imag_dev;
  else
    insist(0);
}

void OmegaWavepacket::copy_weighted_psi_from_host_to_device()
{
  std::cout << " Copy OmegaWavepacket from host to device, Omega: " << omega << std::endl;

  setup_weighted_psi();

  insist(weighted_psi_real && weighted_psi_imag);

  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  if(!weighted_psi_real_dev) 
    checkCudaErrors(hipMalloc(&weighted_psi_real_dev, n1*n2*n_theta*sizeof(double)));
  
  if(!weighted_psi_imag_dev) 
    checkCudaErrors(hipMalloc(&weighted_psi_imag_dev, n1*n2*n_theta*sizeof(double)));
  
  insist(weighted_psi_real_dev && weighted_psi_imag_dev);
  
  checkCudaErrors(hipMemcpyAsync(weighted_psi_real_dev, weighted_psi_real, 
				  n1*n2*n_theta*sizeof(double), 
				  hipMemcpyHostToDevice));
  
  checkCudaErrors(hipMemcpyAsync(weighted_psi_imag_dev, weighted_psi_imag, 
				  n1*n2*n_theta*sizeof(double), 
				  hipMemcpyHostToDevice));
}

void OmegaWavepacket::copy_weighted_psi_from_device_to_host() const
{
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  insist(weighted_psi_real && weighted_psi_imag);
  insist(weighted_psi_real_dev && weighted_psi_imag_dev);
  
  // we can not use unsynchronize version with OpenMP, PCIe can not handle too much data transfer 
  // at the same time, node will be crashed
  
  checkCudaErrors(hipMemcpy(weighted_psi_real, weighted_psi_real_dev, n1*n2*n_theta*sizeof(double), 
			     hipMemcpyDeviceToHost));
  
  checkCudaErrors(hipMemcpy(weighted_psi_imag, weighted_psi_imag_dev, n1*n2*n_theta*sizeof(double), 
			     hipMemcpyDeviceToHost));
}

void OmegaWavepacket::setup_weighted_psi()
{
  if(weighted_psi_real || weighted_psi_imag) {
    insist(weighted_psi_real && weighted_psi_imag);
    return;
  }
  
  const int &omega_min = MatlabData::wavepacket_parameters()->omega_min;
  Vec<RVec> &weighted_wavepackets_real = MatlabData::wavepacket_parameters()->weighted_wavepackets_real;
  Vec<RVec> &weighted_wavepackets_imag = MatlabData::wavepacket_parameters()->weighted_wavepackets_imag;
  
  const int omega_index = omega - omega_min;
  
  weighted_psi_real = weighted_wavepackets_real[omega_index];
  weighted_psi_imag = weighted_wavepackets_imag[omega_index];
}

void OmegaWavepacket::copy_weighted_associated_legendres_from_host_to_device()
{
  if(weighted_associated_legendres_dev) return;

  std::cout << " Copy associated Legendres to device, Omega: " << omega;
  
  const int &n_theta = MatlabData::theta()->n;
  const int &omega_min = MatlabData::wavepacket_parameters()->omega_min;

  const int n_ass_Legs = number_of_associated_Legendres();
  
  const int omega_index = omega - omega_min;
  
  const Vec<RMat> &Legendres = MatlabData::wavepacket_parameters()->weighted_associated_legendres;
  const RMat &ass_Leg = Legendres[omega_index];
  
  insist(ass_Leg.rows() == n_theta && ass_Leg.columns() == n_ass_Legs);
  
  std::cout << ", size: " << n_theta << " " << n_ass_Legs << std::endl;

  checkCudaErrors(hipMalloc(&weighted_associated_legendres_dev, n_theta*n_ass_Legs*sizeof(double)));
  insist(weighted_associated_legendres_dev);
  
  checkCudaErrors(hipMemcpyAsync(weighted_associated_legendres_dev, ass_Leg,
				  n_theta*n_ass_Legs*sizeof(double), hipMemcpyHostToDevice));
}

double OmegaWavepacket::dot_product_with_volume_element(const double *x_dev, const double *y_dev) const
{
  insist(x_dev && y_dev);

  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  const double &dr1 = MatlabData::r1()->dr;
  const double &dr2 = MatlabData::r2()->dr;
  
  double s = 0.0;
  insist(hipblasDdot(cublas_handle, n1*n2*n_theta, x_dev, 1, y_dev, 1, &s) 
	 == HIPBLAS_STATUS_SUCCESS);

  s *= dr1*dr2;
  
  return s;
}

void OmegaWavepacket::calculate_wavepacket_module()
{
  _wavepacket_module_from_real = dot_product_with_volume_element(weighted_psi_real_dev,
								 weighted_psi_real_dev);
  
  _wavepacket_module_from_imag = dot_product_with_volume_element(weighted_psi_imag_dev,
								 weighted_psi_imag_dev);
}

void OmegaWavepacket::calculate_radial_kinetic_add_to_H_weighted_psi_dev() const
{
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;

  double *cufft_tmp_dev = const_cast<double *>(memory_10());
  
  insist(hipfftExecD2Z(cufft_plan_D2Z, (hipfftDoubleReal *) weighted_psi_dev,
		      (hipfftDoubleComplex *) cufft_tmp_dev) == HIPFFT_SUCCESS);
  
  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, (n1/2+1)*n2*n_theta);
  
  _psi_times_kinetic_energy_<<<n_blocks, n_threads, (n1/2+1+n2)*sizeof(double)>>>
    ((Complex *) cufft_tmp_dev, (const Complex *) cufft_tmp_dev, n1, n2, n_theta);
  
  insist(hipfftExecZ2D(cufft_plan_Z2D, (hipfftDoubleComplex *) cufft_tmp_dev,
		      (hipfftDoubleReal *) cufft_tmp_dev) == HIPFFT_SUCCESS);
  
  insist(H_weighted_psi_dev == memory_1());

  n_blocks = cudaUtils::number_of_blocks(n_threads, (n1/2+1)*2*n2*n_theta);
  
  _add_T_radial_weighted_psi_to_H_weighted_psi_
    <<<n_blocks, n_threads>>>(H_weighted_psi_dev, cufft_tmp_dev, n1, n2, n_theta);
}

void OmegaWavepacket::calculate_potential_add_to_H_weighted_psi_dev() const
{
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  
  insist(H_weighted_psi_dev == memory_1());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  const int n_threads = _NTHREADS_;
  const int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n2*n_theta);
  
  _add_potential_weighted_psi_to_H_weighted_psi_<<<n_blocks, n_threads>>>
    (H_weighted_psi_dev, weighted_psi_dev, potential_dev, n1*n2*n_theta);
}

void OmegaWavepacket::setup_work_dev()
{
  if(work_dev) return;
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;

  int size = n1*n2*n_theta;
  
  if(MatlabData::options()->calculate_reaction_probabilities)
    size = std::max(size, 4*n1*n_theta);
  
  std::cout << " Setup local work dev, Omega: " << omega << " " << size << std::endl; 
  
  checkCudaErrors(hipMalloc(&work_dev, size*sizeof(double)));
  insist(work_dev);
}

void OmegaWavepacket::forward_legendre_transform()
{ 
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  insist(weighted_associated_legendres_dev);

  legendre_psi_dev = const_cast<double *>(memory_1());

  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  //const int n_Legs = MatlabData::wavepacket_parameters()->l_max - omega + 1;
  const int n_ass_Legs = number_of_associated_Legendres();

  const double zero = 0.0;
  const double one = 1.0;
  
  double *legendre_psi_dev_ = legendre_psi_dev + omega*n1*n2;

  insist(hipblasDgemm(cublas_handle, HIPBLAS_OP_N, HIPBLAS_OP_N,
		     n1*n2, n_ass_Legs, n_theta, 
		     &one, 
		     weighted_psi_dev, n1*n2,
		     weighted_associated_legendres_dev, n_theta, 
		     &zero,
		     legendre_psi_dev_, n1*n2) == HIPBLAS_STATUS_SUCCESS);
}

void OmegaWavepacket::backward_legendre_transform() const
{ 
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  insist(weighted_associated_legendres_dev);

  insist(legendre_psi_dev == memory_1());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  //const int n_Legs = MatlabData::wavepacket_parameters()->l_max - omega + 1;
  const int n_ass_Legs = number_of_associated_Legendres();
  
  const double zero = 0.0;
  const double one = 1.0;
  
  const double *legendre_psi_dev_ = legendre_psi_dev + omega*n1*n2;
  
  insist(hipblasDgemm(cublas_handle, HIPBLAS_OP_N, HIPBLAS_OP_T,
		     n1*n2, n_theta, n_ass_Legs,
		     &one, 
		     legendre_psi_dev_, n1*n2,
		     weighted_associated_legendres_dev, n_theta, 
		     &zero,
		     weighted_psi_dev, n1*n2) == HIPBLAS_STATUS_SUCCESS);
}

void OmegaWavepacket::T_angle_legendre_psi_to_H_weighted_psi_dev()
{
  insist(T_angle_legendre_psi_dev == memory_10());

  H_weighted_psi_dev = const_cast<double *>(memory_1());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  //const int n_Legs = MatlabData::wavepacket_parameters()->l_max - omega + 1;
  const int n_ass_Legs = number_of_associated_Legendres();
  
  const double zero = 0.0;
  const double one = 1.0;

  const double *T_angle_legendre_psi_dev_ = T_angle_legendre_psi_dev + omega*n1*n2;
  
  insist(hipblasDgemm(cublas_handle, HIPBLAS_OP_N, HIPBLAS_OP_T,
		     n1*n2, n_theta, n_ass_Legs,
		     &one, 
		     T_angle_legendre_psi_dev_, n1*n2,
		     weighted_associated_legendres_dev, n_theta, 
		     &zero,
		     H_weighted_psi_dev, n1*n2) == HIPBLAS_STATUS_SUCCESS);
}

void OmegaWavepacket::copy_T_angle_legendre_psi_to_device_work_dev()
{
  insist(T_angle_legendre_psi_dev == memory_0());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &l_max = MatlabData::wavepacket_parameters()->l_max;
  
  checkCudaErrors(hipMemcpy(device_work_dev, T_angle_legendre_psi_dev, 
			     n1*n2*(l_max+1)*sizeof(double),
			     hipMemcpyDeviceToDevice));
  
  T_angle_legendre_psi_dev = const_cast<double *>(memory_10());
}

void OmegaWavepacket::calculate_T_bend_T_sym_add_to_T_angle_legendre_psi_dev()
{ 
  insist(legendre_psi_dev == memory_1());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  const int &J = MatlabData::wavepacket_parameters()->J;
  
  insist(computation_stream);
  
  T_angle_legendre_psi_dev = const_cast<double *>(memory_0());
  
  checkCudaErrors(hipMemsetAsync(T_angle_legendre_psi_dev, 0, n1*n2*n_theta*sizeof(double),
				  *computation_stream));
  
  //const int n_Legs = MatlabData::wavepacket_parameters()->l_max - omega + 1;
  const int n_ass_Legs = number_of_associated_Legendres();
  
  double *T_angle_legendre_psi_dev_ = T_angle_legendre_psi_dev + omega*n1*n2;
  const double *legendre_psi_dev_ = legendre_psi_dev + omega*n1*n2;
  
  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n2*n_ass_Legs);

  int a = 1; int b = omega;
  
  const int &rot_states = MatlabData::options()->rotational_states;
  if(rot_states == _RotStatesOdd_) {
    a = 2; b = EvolutionUtils::int_to_odd_right(omega);
  } else if(rot_states == _RotStatesEven_) {
    a = 2; b = EvolutionUtils::int_to_even_right(omega);
  } 
  
  _add_T_bend_T_sym_to_T_angle_legendre_psi_dev_<<<n_blocks, n_threads, 
    (n1+n2+1)*sizeof(double), *computation_stream>>>(T_angle_legendre_psi_dev_, 
						     legendre_psi_dev_,
						     n1, n2, n_ass_Legs, J, omega,
						     a, b);
}

void OmegaWavepacket::calculate_T_asym_add_to_T_angle_legendre_psi_dev(const double *psi_dev, 
								       const int omega1) const
{
  insist(omega1 == omega+1 || omega1 == omega-1);
  
  insist(T_angle_legendre_psi_dev == memory_0());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &J = MatlabData::wavepacket_parameters()->J;
  
  const int omega_max = std::max(omega, omega1);
  
  // const int n_Legs = MatlabData::wavepacket_parameters()->l_max - omega_max + 1;
  const int n_ass_Legs = number_of_associated_Legendres();

  double *T_angle_legendre_psi_dev_ = T_angle_legendre_psi_dev + omega_max*n1*n2;
  const double *legendre_psi_dev_ = psi_dev + omega_max*n1*n2;

  int a = 1; int b = omega_max;
  
  const int &rot_states = MatlabData::options()->rotational_states;
  if(rot_states == _RotStatesOdd_) {
    a = 2; b = EvolutionUtils::int_to_odd_right(omega_max);
  } else if(rot_states == _RotStatesEven_) {
    a = 2; b = EvolutionUtils::int_to_even_right(omega_max);
  } 

  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n2*n_ass_Legs);
  
  insist(computation_stream);
  
  _add_T_asym_to_T_angle_legendre_psi_dev_<<<n_blocks, n_threads, n1*sizeof(double),
    *computation_stream>>>(T_angle_legendre_psi_dev_, legendre_psi_dev_,
			   n1, n2, n_ass_Legs, J, omega, omega1, omega_max,
			   a, b);
}


void OmegaWavepacket::calculate_H_weighted_psi_dev()
{
  copy_T_angle_legendre_psi_to_device_work_dev();

  backward_legendre_transform();

  T_angle_legendre_psi_to_H_weighted_psi_dev();

  calculate_radial_kinetic_add_to_H_weighted_psi_dev();
  
  calculate_potential_add_to_H_weighted_psi_dev();
}

void OmegaWavepacket::calculate_energy_and_module()
{   
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  insist(H_weighted_psi_dev == memory_1());
  
  const double module = dot_product_with_volume_element(weighted_psi_dev, weighted_psi_dev);
  const double energy = dot_product_with_volume_element(weighted_psi_dev, H_weighted_psi_dev);
  
  if(weighted_psi_dev == weighted_psi_real_dev) {
    _wavepacket_module_from_real = module;
    _energy_from_real = energy;
  } else {
    _wavepacket_module_from_imag = module;
    _energy_from_imag = energy;
  }
}

void OmegaWavepacket::propagate_with_symplectic_integrator(const int i_step)
{
  insist(weighted_psi_dev == weighted_psi_real_dev || weighted_psi_dev == weighted_psi_imag_dev);
  insist(H_weighted_psi_dev == memory_1());
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  const double &dt = MatlabData::time()->time_step;
  
  //const int &size = SymplecticUtils::coeffients_m6_n4.size;
  const int &size = SymplecticUtils::size();
  
  insist(i_step < size);
  
  if(i_step == size-1) calculate_energy_and_module();

  double coeff = 0.0;
  
  if(weighted_psi_dev == weighted_psi_real_dev) {
    weighted_psi_dev = weighted_psi_imag_dev;
    coeff = -SymplecticUtils::b()[i_step];
  } else {
    weighted_psi_dev = weighted_psi_real_dev;
    coeff = SymplecticUtils::a()[i_step];
  }
  
  coeff *= dt;
  
  insist(hipblasDaxpy(cublas_handle, n1*n2*n_theta, &coeff, 
		     H_weighted_psi_dev, 1,
		     weighted_psi_dev, 1) == HIPBLAS_STATUS_SUCCESS);
  
#if 0
  const double &potential_cutoff = MatlabData::options()->potential_cutoff;
  if(potential_cutoff > _POTENTIAL_CUTOFF_) {
    
    const int &n1 = MatlabData::r1()->n;
    const int &n2 = MatlabData::r2()->n;
    const int &n_theta = MatlabData::theta()->n;
    
    const int n_threads = _NTHREADS_;
    int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n2*n_theta);
    
    _scale_wavepacket_with_potential_cutoff_<<<n_blocks, n_threads>>>(weighted_psi_dev,
								      potential_dev,
								      potential_cutoff,
								      n1*n2*n_theta);
  }
#endif
}

void OmegaWavepacket::dump_wavepacket() const
{
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n2*n_theta);

  _dump_wavepacket_<<<n_blocks, n_threads>>>(weighted_psi_real_dev, n1, n2, n_theta);
  _dump_wavepacket_<<<n_blocks, n_threads>>>(weighted_psi_imag_dev, n1, n2, n_theta);
}

void OmegaWavepacket::setup_reaction_probabilities()
{
  if(!MatlabData::options()->calculate_reaction_probabilities) return;
  
  if(crp) return;
  
  std::cout << " Setup reaction probabilities, Omega: " << omega << std::endl;
  
  crp = new ReactionProbabilities(this);
  insist(crp);
}

void OmegaWavepacket::calculate_reaction_probabilities(const int calculate)
{ if(crp) crp->calculate_reaction_probabilities(calculate); }

const double *OmegaWavepacket::reaction_probabilities() const
{ return crp ? (const double *) crp->reaction_probabilities : 0; }



