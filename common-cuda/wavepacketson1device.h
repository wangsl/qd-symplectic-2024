
#ifndef WAVEPACKETS_ON_SINGLE_DEVICE
#define WAVEPACKETS_ON_SINGLE_DEVICE

#include <cublas_v2.h>
#include <cufft.h>

#include "omegawavepacket.h"
#include "vecbase.h"
#include "rmat.h"

class WavepacketsOnSingleDevice
{
public:
  WavepacketsOnSingleDevice(const int device_index,
			    const int omega_start,
			    const int n_omegas);

  ~WavepacketsOnSingleDevice() { destroy_data_on_device(); }

  void print();

  void setup_neighbours(const WavepacketsOnSingleDevice *left, 
			const WavepacketsOnSingleDevice *right);

  void setup_device_work_dev_and_copy_streams_events();

  void propagate_with_symplectic_integrator(const int istep);

  void copy_weighted_psi_from_device_to_host();

  void dump_wavepackets() const;

  double module() const { return _module; }
  double total_energy() const { return _total_energy; }

  void calculate_reaction_probabilities(const int calculate);

  RVec reaction_probabilities;

  void test();

private:

  void forward_legendre_transform_and_copy_data_to_neighbour_devices(const int part);

  void calculate_T_asym_add_to_T_angle_legendre_psi_dev();

  void calculate_H_weighted_psi_dev(const int part);

private:
  
  Vec<OmegaWavepacket *> omega_wavepackets;

  int _device_index;
  int omega_start;
  int n_omegas;

  double *potential_dev;

  double *device_work_dev;

  cublasHandle_t cublas_handle;
  int _has_created_cublas_handle;
  
  cufftHandle cufft_plan_D2Z;
  cufftHandle cufft_plan_Z2D;
  int _has_cufft_plans;
  
  const WavepacketsOnSingleDevice *left;
  const WavepacketsOnSingleDevice *right;
  
  double *omega_wavepacket_from_left_device;
  double *omega_wavepacket_from_right_device;

  cudaStream_t *computation_stream;
  cudaStream_t *data_copy_stream;
  
  cudaEvent_t *computation_event_for_left;
  cudaEvent_t *computation_event_for_right;
  
  cudaEvent_t *copy_to_left_event;
  cudaEvent_t *copy_to_right_event;

  double _module;
  double _total_energy;
  
  int device_index() const { return _device_index; }
  int current_device_index() const;
  void setup_device() const;

  void setup_data_on_device();
  void destroy_data_on_device();

  void setup_constant_memory_on_device();

  void setup_cublas_handle();
  void destroy_cublas_handle();

  void setup_cufft_plans();
  void destroy_cufft_plans();

  void setup_potential_on_device();
  void setup_omega_wavepackets();

  void setup_computation_stream_and_event();
  void destroy_streams_and_events();

  int copy_to_left_event_query() const;
  int copy_to_right_event_query() const;
  int ready_to_receive_data() const;

  void copy_numerical_gradient_coefficients_to_device() const;
  void copy_reaction_probabity_energies_to_device() const;
};

#endif /* WAVEPACKETS_ON_SINGLE_DEVICE */
