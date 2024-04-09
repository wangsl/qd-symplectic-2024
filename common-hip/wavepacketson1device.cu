
#include <iostream>

#include "wavepacketson1device.h"
#include "cudaUtils.h"
#include "matlabUtils.h"
#include "matlabData.h"
#include "evolutionUtils.h"

#include "evolutionAux.cu"

/***
 * https://github.com/mohamso/icpads14
 * https://raw.githubusercontent.com/mohamso/icpads14/master/4/omp/src/Async.c
 * DOI: 10.1109/PADSW.2014.7097919
 ***/

__constant__ EvolutionUtils::RadialCoordinate r1_dev;
__constant__ EvolutionUtils::RadialCoordinate r2_dev;
__constant__ double energies_dev[_EnergiesMaxSize_];
__constant__ double potential_cutoff;
//__constant__ double kinetic_cutoff;

WavepacketsOnSingleDevice::
WavepacketsOnSingleDevice(const int device_index_,
			  const int omega_start_,
			  const int n_omegas_) :
  _device_index(device_index_),
  omega_start(omega_start_),
  n_omegas(n_omegas_),
  potential_dev(0),
  device_work_dev(0),
  omega_wavepacket_from_left_device(0),
  omega_wavepacket_from_right_device(0),
  _has_created_cublas_handle(0),
  _has_cufft_plans(0),
  computation_stream(0), 
  computation_event_for_left(0),
  computation_event_for_right(0),
  copy_to_left_event(0),
  copy_to_right_event(0),
  data_copy_stream(0),
  left(0), right(0)
{ 
  insist(_device_index >= 0);
  setup_data_on_device();
}

int WavepacketsOnSingleDevice::current_device_index() const
{
  int dev_index = -1;
  checkCudaErrors(hipGetDevice(&dev_index));
  return dev_index;
}

void WavepacketsOnSingleDevice::setup_device() const
{
  if(current_device_index() != device_index()) 
    checkCudaErrors(hipSetDevice(device_index()));
}

void WavepacketsOnSingleDevice::setup_data_on_device()
{
  setup_device();

  std::cout << " Setup data on device: " << device_index() << std::endl;

  setup_constant_memory_on_device();
  
  setup_computation_stream_and_event();

  setup_cublas_handle();
  setup_cufft_plans();

  setup_potential_on_device();
  setup_omega_wavepackets();
}

void WavepacketsOnSingleDevice::destroy_data_on_device()
{ 
  setup_device();
  
  std::cout << " Destroy data on device: " << device_index() << std::endl;

  for(int i = 0; i < omega_wavepackets.size(); i++) 
    if(omega_wavepackets[i]) { delete omega_wavepackets[i]; omega_wavepackets[i] = 0; }
  omega_wavepackets.resize(0);

  _CUDA_FREE_(potential_dev);
  _CUDA_FREE_(device_work_dev);
  
  destroy_cublas_handle();
  destroy_cufft_plans();

  destroy_streams_and_events();

  reaction_probabilities.resize(0);

  left = 0;
  right = 0;
}

void WavepacketsOnSingleDevice::setup_potential_on_device()
{
  if(potential_dev) return;

  std::cout << " Allocate and copy potential on device: " << current_device_index() << std::endl;
  
  const double *potential = MatlabData::potential();
  insist(potential);
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  
  checkCudaErrors(hipMalloc(&potential_dev, n1*n2*n_theta*sizeof(double)));
  insist(potential_dev);
  
  checkCudaErrors(hipMemcpyAsync(potential_dev, potential, n1*n2*n_theta*sizeof(double),
				  hipMemcpyHostToDevice));
}

void WavepacketsOnSingleDevice::setup_cublas_handle()
{
  if(_has_created_cublas_handle) return;

  std::cout << " Setup cuBLAS handle on device: " << current_device_index() << std::endl;
  
  insist(hipblasCreate(&cublas_handle) == HIPBLAS_STATUS_SUCCESS);
  
  _has_created_cublas_handle = 1;
}

void WavepacketsOnSingleDevice::destroy_cublas_handle()
{ 
  if(!_has_created_cublas_handle) return;
  
  std::cout << " Destroy cuBLAS handle on device: " << current_device_index() << std::endl;

  insist(hipblasDestroy(cublas_handle) == HIPBLAS_STATUS_SUCCESS);

  _has_created_cublas_handle = 0;
}

void WavepacketsOnSingleDevice::setup_cufft_plans()
{
  if(_has_cufft_plans) return;

  std::cout << " Setup cuFFT handles on device: " << current_device_index() << std::endl;

  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;

  /* wavepacket psi is from Matlab in column-major format, 
   * while cuFFT is using row-major format,
   * so to switch dimensions, after D2Z FFT, the output data is { n2, n1/2+1 }, 
   * it is still in column-major format
   */
  const int dims [] = { n2, n1 };
  
  insist(hipfftPlanMany(&cufft_plan_D2Z, 2, const_cast<int *>(dims), 
		       NULL, 1, n1*n2,
		       NULL, 1, n1*n2,
		       HIPFFT_D2Z, n_theta) == HIPFFT_SUCCESS);

  cudaUtils::cufft_work_size(cufft_plan_D2Z, "D2Z");

  insist(hipfftPlanMany(&cufft_plan_Z2D, 2, const_cast<int *>(dims), 
		       NULL, 1, n1*n2,
		       NULL, 1, n1*n2,
		       HIPFFT_Z2D, n_theta) == HIPFFT_SUCCESS);

  cudaUtils::cufft_work_size(cufft_plan_Z2D, "Z2D");
  
  _has_cufft_plans = 1;
}

void WavepacketsOnSingleDevice::destroy_cufft_plans()
{ 
  if(!_has_cufft_plans) return;

  std::cout << " Destroy cuFFT handles on device: " << current_device_index() << std::endl;

  insist(hipfftDestroy(cufft_plan_D2Z) == HIPFFT_SUCCESS);
  insist(hipfftDestroy(cufft_plan_Z2D) == HIPFFT_SUCCESS);

  _has_cufft_plans = 0;
}

void WavepacketsOnSingleDevice::setup_omega_wavepackets()
{
  insist(omega_wavepackets.size() == 0);
  
  omega_wavepackets.resize(n_omegas, 0);
  
  for(int i = 0; i < n_omegas; i++) {
    omega_wavepackets[i] = new OmegaWavepacket(i+omega_start, potential_dev, 
					       cublas_handle, cufft_plan_D2Z, cufft_plan_Z2D,
					       computation_stream,
					       device_work_dev);
    insist(omega_wavepackets[i]);
  }
}

void WavepacketsOnSingleDevice::setup_constant_memory_on_device()
{
  std::cout << " Setup constant memory on device: " << current_device_index() << std::endl;

  EvolutionUtils::copy_radial_coordinate_to_device(r1_dev, MatlabData::r1());
  EvolutionUtils::copy_radial_coordinate_to_device(r2_dev, MatlabData::r2());

  checkCudaErrors(hipMemcpyToSymbolAsync(HIP_SYMBOL(potential_cutoff),
					  &MatlabData::options()->potential_cutoff,
					  sizeof(double),
					  0,
					  hipMemcpyHostToDevice,
					  0
					  ));

#if 0
  checkCudaErrors(hipMemcpyToSymbolAsync(HIP_SYMBOL(kinetic_cutoff),
					  &MatlabData::options()->kinetic_cutoff,
					  sizeof(double)));
#endif

  copy_numerical_gradient_coefficients_to_device();
  copy_reaction_probabity_energies_to_device();
}

void WavepacketsOnSingleDevice::setup_device_work_dev_and_copy_streams_events()
{
  setup_device();

  if(device_work_dev) return;
  
  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int &n_theta = MatlabData::theta()->n;
  const int &l_max = MatlabData::wavepacket_parameters()->l_max;
  
  long size = 0;

  if(left) size += n1*n2*(l_max+1);
  if(right) size += n1*n2*(l_max+1);
  
  size = std::max(size, (n1/2+1)*2L*n2*n_theta);
  
  std::cout << " Setup device work on device: " << current_device_index() 
	    << " " << size << " " << size*sizeof(double)/1024.0/1024.0 << std::endl;
  
  checkCudaErrors(hipMalloc(&device_work_dev, size*sizeof(double)));
  insist(device_work_dev);
  
  long current = 0;
  
  if(left || right) _CUDA_STREAM_CREATE_(data_copy_stream);
  
  if(left && !omega_wavepacket_from_left_device) {
    omega_wavepacket_from_left_device = device_work_dev + current;
    current += n1*n2*(l_max+1);
    
    std::cout << " Setup wavepacket from left on device: " << current_device_index() 
	      << " " << omega_wavepacket_from_left_device << std::endl;

    _CUDA_EVENT_CREATE_(computation_event_for_left);
    _CUDA_EVENT_CREATE_(copy_to_left_event);
  }
  
  if(right && !omega_wavepacket_from_right_device) {
    omega_wavepacket_from_right_device = device_work_dev + current;
    current += n1*n2*(l_max+1);
    
    std::cout << " Setup wavepacket from right on device: " << current_device_index() 
	      << " " << omega_wavepacket_from_right_device << std::endl;

    _CUDA_EVENT_CREATE_(computation_event_for_right);
    _CUDA_EVENT_CREATE_(copy_to_right_event);
  }
}

void WavepacketsOnSingleDevice::setup_computation_stream_and_event()
{ 
  if(computation_stream) return;
  
  std::cout << " Setup computation stream on device: " << current_device_index() << std::endl;
  
  _CUDA_STREAM_CREATE_(computation_stream);
}

void WavepacketsOnSingleDevice::destroy_streams_and_events()
{ 
  if(cudaUtils::n_devices() == 1) return;

  std::cout << " Destroy streams and events on device: " << device_index() << std::endl;
  
  _CUDA_STREAM_DESTROY_(computation_stream);
  _CUDA_STREAM_DESTROY_(data_copy_stream);

  _CUDA_EVENT_DESTROY_(computation_event_for_left);
  _CUDA_EVENT_DESTROY_(computation_event_for_right);
  _CUDA_EVENT_DESTROY_(copy_to_left_event);
  _CUDA_EVENT_DESTROY_(copy_to_right_event);
}

void WavepacketsOnSingleDevice::setup_neighbours(const WavepacketsOnSingleDevice *left_, 
						 const WavepacketsOnSingleDevice *right_)
{
  setup_device();
  left = left_; 
  right = right_;
  std::cout << " Neighbours on device: " << current_device_index()
	    << ", pointers: " << this << " " << left << " " << right << std::endl;
}

void WavepacketsOnSingleDevice::
forward_legendre_transform_and_copy_data_to_neighbour_devices(const int part)
{ 
  insist(part == _RealPart_ || part == _ImagPart_);

  setup_device();

  const int &n1 = MatlabData::r1()->n;
  const int &n2 = MatlabData::r2()->n;
  const int n_Legs = MatlabData::wavepacket_parameters()->l_max + 1;
  
  insist(computation_stream);
  insist(hipblasSetStream(cublas_handle, *computation_stream) == HIPBLAS_STATUS_SUCCESS);
  
  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->setup_weighted_psi_dev(part);
  
  OmegaWavepacket *wp_for_left = 0;
  OmegaWavepacket *wp_for_right = 0;
  
  int wp_start = 0;
  int wp_end = n_omegas;

  if(n_omegas == 1) {
    
    if(left || right) 
      omega_wavepackets[0]->forward_legendre_transform();
    
    if(left) {
      insist(computation_event_for_left);
      checkCudaErrors(hipEventRecord(*computation_event_for_left, *computation_stream));
      wp_for_left = omega_wavepackets[0];
      wp_start = 1;
    }
    
    if(right) {
      insist(computation_event_for_right);
      checkCudaErrors(hipEventRecord(*computation_event_for_right, *computation_stream));
      wp_for_right = omega_wavepackets[0];
      wp_end = n_omegas - 1;
    }

  } else {
    
    if(left) {
      insist(computation_event_for_left);
      wp_for_left = omega_wavepackets[0];
      wp_for_left->forward_legendre_transform();
      checkCudaErrors(hipEventRecord(*computation_event_for_left, *computation_stream));
      wp_start = 1;
    }

    if(right) {
      insist(computation_event_for_right);
      wp_for_right = omega_wavepackets[n_omegas-1];
      wp_for_right->forward_legendre_transform();
      checkCudaErrors(hipEventRecord(*computation_event_for_right, *computation_stream));
      wp_end = n_omegas - 1;
    }
  }
  
  int sent_to_left = left ? 0 : 1;
  int sent_to_right = right ? 0 : 1;
  
  while(true) {
    
    if(left && !sent_to_left) {
      
      if(left->ready_to_receive_data()) {
	
	checkCudaErrors(hipStreamWaitEvent(*data_copy_stream, *computation_event_for_left, 0));
	  
	checkCudaErrors(hipMemcpyPeerAsync(left->omega_wavepacket_from_right_device, left->device_index(),
					    wp_for_left->legendre_psi_dev_(), device_index(),
					    n1*n2*n_Legs*sizeof(double), *data_copy_stream));
	
	checkCudaErrors(hipEventRecord(*copy_to_left_event, *data_copy_stream));
	
	sent_to_left = 1;
      }
    }
    
    if(right && !sent_to_right) {
      
      if(right->ready_to_receive_data()) {
	
	checkCudaErrors(hipStreamWaitEvent(*data_copy_stream, *computation_event_for_right, 0));
	
	checkCudaErrors(hipMemcpyPeerAsync(right->omega_wavepacket_from_left_device, right->device_index(),
					    wp_for_right->legendre_psi_dev_(), device_index(),
					    n1*n2*n_Legs*sizeof(double), *data_copy_stream));
	
	checkCudaErrors(hipEventRecord(*copy_to_right_event, *data_copy_stream));
	
	sent_to_right = 1;
      }
    }
    
    if(sent_to_left && sent_to_right) break;
  }
  
  for(int i = wp_start; i < wp_end; i++) 
    omega_wavepackets[i]->forward_legendre_transform();
}

void WavepacketsOnSingleDevice::calculate_T_asym_add_to_T_angle_legendre_psi_dev()
{
  setup_device();

  const OmegaWavepacket *wp = 0;

  for(int i = 0; i < n_omegas; i++) {
    if(i > 0) {
      wp = omega_wavepackets[i-1];
      omega_wavepackets[i]->calculate_T_asym_add_to_T_angle_legendre_psi_dev(wp->legendre_psi_dev_(),
									     wp->omega_());
    }
    
    if(i < n_omegas-1) {
      wp = omega_wavepackets[i+1];
      omega_wavepackets[i]->calculate_T_asym_add_to_T_angle_legendre_psi_dev(wp->legendre_psi_dev_(),
									     wp->omega_());
    }
  }

#pragma omp barrier
  
  if(right) {
    insist(right->copy_to_left_event);

    checkCudaErrors(hipStreamWaitEvent(*computation_stream, *right->copy_to_left_event, 0));
    
    wp = omega_wavepackets[n_omegas-1];
    wp->calculate_T_asym_add_to_T_angle_legendre_psi_dev(omega_wavepacket_from_right_device,
							 wp->omega_()+1);
  }

  if(left) {
    insist(left->copy_to_right_event);
    
    checkCudaErrors(hipStreamWaitEvent(*computation_stream, *left->copy_to_right_event, 0));
    
    wp = omega_wavepackets[0];
    wp->calculate_T_asym_add_to_T_angle_legendre_psi_dev(omega_wavepacket_from_left_device,
							 wp->omega_()-1);
  }
}

void WavepacketsOnSingleDevice::calculate_H_weighted_psi_dev(const int part)
{
  forward_legendre_transform_and_copy_data_to_neighbour_devices(part);

  for(int i = 0; i < n_omegas; i++)
    omega_wavepackets[i]->calculate_T_bend_T_sym_add_to_T_angle_legendre_psi_dev();

  calculate_T_asym_add_to_T_angle_legendre_psi_dev();

  checkCudaErrors(hipDeviceSynchronize());
  insist(hipblasSetStream(cublas_handle, NULL) == HIPBLAS_STATUS_SUCCESS);
  
  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->calculate_H_weighted_psi_dev();
}

void WavepacketsOnSingleDevice::propagate_with_symplectic_integrator(const int i_step)
{
  setup_device();

  calculate_H_weighted_psi_dev(_RealPart_);

  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->propagate_with_symplectic_integrator(i_step);
  
  checkCudaErrors(hipDeviceSynchronize());
  
#pragma omp barrier
  
  calculate_H_weighted_psi_dev(_ImagPart_);
  
  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->propagate_with_symplectic_integrator(i_step);
}

void WavepacketsOnSingleDevice::print()
{
  setup_device();
  _module = 0.0;
  _total_energy = 0.0;
  for(int i = 0; i < n_omegas; i++) {
    _module += omega_wavepackets[i]->wavepacket_module();
    _total_energy += omega_wavepackets[i]->energy();
    
    std::cout << " " << omega_wavepackets[i]->omega_()
	      << " " << omega_wavepackets[i]->wavepacket_module()
	      << " " << omega_wavepackets[i]->energy()
	      << std::endl;
  }
}

void WavepacketsOnSingleDevice::copy_weighted_psi_from_device_to_host()
{
  setup_device();
  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->copy_weighted_psi_from_device_to_host();
}

void WavepacketsOnSingleDevice::dump_wavepackets() const
{
  setup_device();
  for(int i = 0; i < n_omegas; i++) 
    omega_wavepackets[i]->dump_wavepacket();
}

int WavepacketsOnSingleDevice::copy_to_left_event_query() const
{
  if(copy_to_left_event) 
    return hipEventQuery(*copy_to_left_event) == hipSuccess ? 1 : 0;
  else 
    return 1;
}

int WavepacketsOnSingleDevice::copy_to_right_event_query() const
{
  if(copy_to_right_event) 
    return hipEventQuery(*copy_to_right_event) == hipSuccess ? 1 : 0;
  else 
    return 1;
}

int WavepacketsOnSingleDevice::ready_to_receive_data() const
{
  int left_ok = 1;
  int right_ok = 1;
  
  if(left) left_ok = left->copy_to_right_event_query();
  if(right) right_ok = right->copy_to_left_event_query();

  return left_ok && right_ok ? 1 : 0;
}

void WavepacketsOnSingleDevice::copy_numerical_gradient_coefficients_to_device() const
{
  if(!MatlabData::options()->calculate_reaction_probabilities) return;

  std::cout << " Copy numerical gradient coefficients on device: " 
	    << current_device_index() << std::endl;

  const int &n_points = MatlabData::crp_parameters()->n_gradient_points;
  
  Num1ststGradient::copy_gradient_coefficients_to_device(n_points);
}

void WavepacketsOnSingleDevice::test()
{ 
  setup_device();

  const int &n_points = MatlabData::crp_parameters()->n_gradient_points;

  std::cout << " Print " << n_points  << " gradient coefficients on device: " 
	    << current_device_index() << std::endl;
  
  _print_gradient_coeffients_<<<1,1>>>(n_points/2);

  const int &n_energies = MatlabData::crp_parameters()->n_energies;
  
  _print_energies_<<<1,1>>>(n_energies);

  checkCudaErrors(hipDeviceSynchronize());
}

void WavepacketsOnSingleDevice::copy_reaction_probabity_energies_to_device() const
{
  if(!MatlabData::options()->calculate_reaction_probabilities) return;
  
  std::cout << " Copy reaction probability enegies to device: " 
	    << current_device_index() << std::endl;
  
  const int &n_energies = MatlabData::crp_parameters()->n_energies;
  const RVec &energies = MatlabData::crp_parameters()->energies;
  
  size_t size = 0;
  checkCudaErrors(hipGetSymbolSize(&size, HIP_SYMBOL(energies_dev)));
  insist(size > n_energies*sizeof(double));
  
  checkCudaErrors(hipMemcpyToSymbolAsync(HIP_SYMBOL(energies_dev), energies,
					  n_energies*sizeof(double), 0, 
					  hipMemcpyHostToDevice));
}

void WavepacketsOnSingleDevice::calculate_reaction_probabilities(const int calculate)
{
  if(!MatlabData::options()->calculate_reaction_probabilities) return;

  setup_device();

  for(int i = 0; i < n_omegas; i++)
    omega_wavepackets[i]->calculate_reaction_probabilities(calculate);

  if(calculate) {
    const int &n_energies = MatlabData::crp_parameters()->n_energies;
    
    reaction_probabilities.resize(n_energies, 0);
    reaction_probabilities.zeros();
    
    for(int i = 0; i < n_omegas; i++) {
      const RVec crp(n_energies, const_cast<double *>(omega_wavepackets[i]->reaction_probabilities()));
      reaction_probabilities += crp;
    }
  }
}
