
#ifndef REACTION_PROBABILITIES_H
#define REACTION_PROBABILITIES_H

#include "complex.h"
#include "omegawavepacket.h"

class ReactionProbabilities
{
  
public:
  
  ReactionProbabilities(const OmegaWavepacket *wavepacket);

  ~ReactionProbabilities();

  RVec reaction_probabilities;

  void calculate_reaction_probabilities(const int calculate);
  
private:
  
  const OmegaWavepacket *wavepacket;
  Complex *fai_on_surface_dev;
  Complex *d_fai_on_surface_dev;
  
  double *psi_real_dev;
  double *d_psi_real_dev;
  double *psi_imag_dev;
  double *d_psi_imag_dev;
  
  void setup_data_on_device();

  void calculate_psi_gradients_on_dividing_surface();
  void psi_time_to_fai_energy_on_dividing_surface() const;
  void calculate_reaction_probabilities();
};

#endif /* REACTION_PROBABILITIES_H */
