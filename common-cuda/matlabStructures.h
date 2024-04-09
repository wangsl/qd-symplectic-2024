
/* $Id$ */

#ifndef MATLAB_STRUCTURES_H
#define MATLAB_STRUCTURES_H

#include <iostream>
using namespace std;

#include <mex.h>
#include "rmat.h"
#include "matlabUtils.h"

class RadialCoordinate
{
public:
  
  const int &n; // out
  const double &left; // out
  const double &dr; // out
  const double &mass; // out
  RVec dump;
  
  RadialCoordinate(const mxArray *mx);

  ~RadialCoordinate() { if(mx) mx = 0; }
  
private:

  const mxArray *mx;
  
  // to prevent assigment and copy operation
  RadialCoordinate(const RadialCoordinate &);
  RadialCoordinate & operator =(const RadialCoordinate &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const RadialCoordinate &c);
  void write_fields(ostream &s) const;
};

class AngleCoordinate
{
public:
  
  const int &n; // out
  RVec x; // out
  RVec w; // out
  
  AngleCoordinate(const mxArray *mx);

  ~AngleCoordinate() { if(mx) mx = 0; }
  
private:

  const mxArray *mx;
  
  // to prevent assigment and copy operation
  AngleCoordinate(const AngleCoordinate &);
  AngleCoordinate & operator =(const AngleCoordinate &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const AngleCoordinate &c);
  void write_fields(ostream &s) const;
};

class EvolutionTime
{
public:

  const int &total_steps; // out
  const double &time_step; // out
  int &steps; // out
  
  EvolutionTime(const mxArray *mx);
  
  ~EvolutionTime() { if(mx) mx = 0; }
  
private:
  
  const mxArray *mx;

  EvolutionTime(const EvolutionTime &);
  EvolutionTime & operator =(const EvolutionTime &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const EvolutionTime &c);
  void write_fields(ostream &s) const;
};

class SICoefficients
{
public:
  const int &m;
  RVec a; // out
  RVec b; // out

  SICoefficients(const mxArray *mx);
  ~SICoefficients() { if(mx) mx = 0; }

public:

  const mxArray *mx;
  
  // to prevent assigment and copy operation
  SICoefficients(const SICoefficients &);
  SICoefficients & operator =(const SICoefficients &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const SICoefficients &c);
  void write_fields(ostream &s) const;
};

class Options
{
public:
  
  char *wave_to_matlab; // out
  const int &steps_to_copy_psi_from_device_to_host; // out
  const double &potential_cutoff; // out
  const int &calculate_reaction_probabilities; // out
  int &rotational_states; // out
  const double &converged_criterion_with_module; // out
  const double &catastrophe_criterion_with_module; // out
  
  Options(const mxArray *mx);
  ~Options();
  
private:
  
  const mxArray *mx;

  Options(const Options &);
  Options & operator =(const Options &);

  friend ostream & operator <<(ostream &s, const Options &c);
  void write_fields(ostream &s) const;
};

class WavepacketParameters
{
public:

  WavepacketParameters(const mxArray *mx);

  ~WavepacketParameters() { if(mx) mx = 0; }  

  const int &J; // out
  const int &parity; // out
  const int &l_max; // out
  const int &omega_min; // out
  const int &omega_max; // out
  
  Vec<RMat> weighted_associated_legendres;
  Vec<RVec> weighted_wavepackets_real;
  Vec<RVec> weighted_wavepackets_imag;

  int legendre_wavepackets_size() { return l_max - omega_min + 1; }

private:

  const mxArray *mx;

  void setup_weighted_associated_legendres();
  void setup_weighted_wavepackets();

  WavepacketParameters(const WavepacketParameters &);
  WavepacketParameters & operator =(const WavepacketParameters &);

  friend ostream & operator <<(ostream &s, const WavepacketParameters &c);
  void write_fields(ostream &s) const;
};

class CRPParameters
{
public:

  CRPParameters(const mxArray *mx);

  ~CRPParameters() { if(mx) mx = 0; }

  const int &n_dividing_surface; // out
  const int &n_gradient_points; // out
  const int &n_energies; // out
  
  RVec energies; 
  RVec eta_sq; 
  RVec CRP; 

private:
  const mxArray *mx;
  
  CRPParameters(const CRPParameters &);
  CRPParameters & operator =(const CRPParameters &);
  
  friend ostream & operator <<(ostream &s, const CRPParameters &c);
  void write_fields(ostream &s) const;
};

#endif /* MATLAB_STRUCTURES_H */
