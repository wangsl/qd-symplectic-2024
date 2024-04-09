
#ifndef CUDA_OPENMP_MD_H
#define CUDA_OPENMP_MD_H

#include "wavepacketson1device.h"
#include "vecbase.h"

class CUDAOpenmpMD
{
public:

  CUDAOpenmpMD();

  ~CUDAOpenmpMD();

  void time_evolution();
  
  void test();

private:
  
  int _n_devices;
  
  Vec<WavepacketsOnSingleDevice *> wavepackets_on_single_device;

  int n_devices() const { return _n_devices; }

  void setup_n_devices();
  void devices_synchoronize();
  void devices_memory_usage() const;
  void reset_devices();
  
  void setup_wavepackets_on_single_device();
  void destroy_wavepackets_on_single_device();

  void enable_peer_to_peer_access() const;
  void disable_peer_to_peer_access() const;

  void setup_devices_neighbours() const;

  void setup_device_work_dev_on_devices() const;

  void copy_weighted_psi_from_device_to_host();

  void dump_wavepackets() const;

  void calculate_reaction_probabilities(const int calculate);
};

#endif /* CUDA_OPENMP_MD_H */
