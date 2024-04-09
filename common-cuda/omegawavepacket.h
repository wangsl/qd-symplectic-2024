
#ifndef OMEGA_WAVEPACKET_H
#define OMEGA_WAVEPACKET_H

#include <cublas_v2.h>
#include <cufft.h>

class ReactionProbabilities;

class OmegaWavepacket
{
  friend ReactionProbabilities;

public:

  OmegaWavepacket(int omega,
		  const double *potential_dev, 
		  cublasHandle_t &cublas_handle,
		  cufftHandle &cufft_plan_D2Z,
		  cufftHandle &cufft_plan_Z2D,
		  cudaStream_t * &computation_stream,
		  double * &device_work_dev
		  );

  ~OmegaWavepacket();

  void calculate_wavepacket_module();

  void calculate_T_bend_T_sym_add_to_T_angle_legendre_psi_dev();

  void calculate_T_asym_add_to_T_angle_legendre_psi_dev(const double *psi_dev, 
							const int omega1) const;
  
  void calculate_H_weighted_psi_dev();

  void propagate_with_symplectic_integrator(const int istep);

  void dump_wavepacket() const;
  
  double wavepacket_module() const 
  { return _wavepacket_module_from_real + _wavepacket_module_from_imag; }

  double energy() const
  { return _energy_from_real + _energy_from_imag; }

  double wavepacket_module_legendre() const 
  { return _wavepacket_module_from_real_legendre + _wavepacket_module_from_imag_legendre; }

  double kinetic_energy() const
  { return _kinetic_energy_from_real + _kinetic_energy_from_imag; }

  double potential_energy() const
  { return _potential_energy_from_real + _potential_energy_from_imag; }

  void setup_weighted_psi_dev(const int part);

  void forward_legendre_transform();

  void copy_weighted_psi_from_device_to_host() const;

  const double *legendre_psi_dev_() const { return work_dev; }
  const int &omega_() const { return omega; }

  void calculate_reaction_probabilities(const int calculate);

  const double *reaction_probabilities() const;

private:

  const int omega;
  
  double *weighted_psi_real;
  double *weighted_psi_imag;
  
  const double *potential_dev;
  double *work_dev;
  
  double *weighted_psi_real_dev;
  double *weighted_psi_imag_dev;
  double *weighted_associated_legendres_dev;
  
  double *weighted_psi_dev;
  double *legendre_psi_dev;

  double *T_angle_legendre_psi_dev;
  double *H_weighted_psi_dev;

  double * &device_work_dev;

  cublasHandle_t &cublas_handle;
  cufftHandle &cufft_plan_D2Z;
  cufftHandle &cufft_plan_Z2D;
  cudaStream_t * &computation_stream;

  double _wavepacket_module_from_real;
  double _wavepacket_module_from_imag;

  double _energy_from_real;
  double _energy_from_imag;

  double _wavepacket_module_from_real_legendre;
  double _wavepacket_module_from_imag_legendre;
  
  double _kinetic_energy_from_real;
  double _kinetic_energy_from_imag;

  double _potential_energy_from_real;
  double _potential_energy_from_imag;

  ReactionProbabilities *crp;

  void setup_weighted_psi();
  void copy_weighted_psi_from_host_to_device();
  
  void copy_weighted_associated_legendres_from_host_to_device();
  
  void setup_work_dev();

  void backward_legendre_transform() const;

  void calculate_radial_kinetic_add_to_H_weighted_psi_dev() const;

  void calculate_radial_kinetic_add_to_H_weighted_psi_dev_2() const;

  void calculate_potential_add_to_H_weighted_psi_dev() const;

  void copy_T_angle_legendre_psi_to_device_work_dev();
  
  void T_angle_legendre_psi_to_H_weighted_psi_dev();

  void calculate_energy_and_module();

  double dot_product_with_volume_element(const double *x_dev, const double *y_dev) const;
  //double dot_product_with_volume_element_for_legendres(const double *x_dev, const double *y_dev) const;

  const double *memory_0() const { return weighted_psi_dev; }

  // fixed memory address
  const double *memory_1() const { return work_dev; }
  const double *memory_10() const { return device_work_dev; }

  void setup_reaction_probabilities();

  int number_of_associated_Legendres() const;
  int minimum_Legendres_order() const;
  int l_minimum() const { return minimum_Legendres_order(); }
};

#endif /* OMEGA_WAVEPACKET_H */
