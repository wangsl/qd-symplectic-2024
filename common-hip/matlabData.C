
#include "matlabData.h"

static const RadialCoordinate *_r1 = 0;
static const RadialCoordinate *_r2 = 0;
static const AngleCoordinate *_theta = 0;
static const double *_potential = 0;
static EvolutionTime *_time = 0;
static const Options *_options = 0;
static WavepacketParameters *_wavepacket_parameters = 0;
static CRPParameters *_crp_parameters = 0;
static const SICoefficients *_si_coefficients = 0;

// SI coefficients
const SICoefficients *MatlabData::si_coefficients() { return _si_coefficients; }
void MatlabData::si_coefficients(const SICoefficients *si)
{ insist(si && !_si_coefficients); _si_coefficients = si; }

// r1
const RadialCoordinate *MatlabData::r1() { return _r1; }
void MatlabData::r1(const RadialCoordinate *r) { insist(r && !_r1); _r1 = r; }

// r2
const RadialCoordinate *MatlabData::r2() { return _r2; }
void MatlabData::r2(const RadialCoordinate *r) { insist(r && !_r2); _r2 = r; }

// theta
const AngleCoordinate *MatlabData::theta() { return _theta; }
void MatlabData::theta(const AngleCoordinate *th) { insist(th && !_theta); _theta = th; }

// potential
const double *MatlabData::potential() { return _potential; }
void MatlabData::potential(const double *p) { insist(p && !_potential); _potential = p; }

// evolution time
EvolutionTime *MatlabData::time() { return _time; }
void MatlabData::time(EvolutionTime *t) { insist(t && !_time); _time = t; }

// options
const Options *MatlabData::options() { return _options; }
void MatlabData::options(const Options *opt) { insist(opt && !_options); _options = opt; }

// wavepacket parameters
WavepacketParameters *MatlabData::wavepacket_parameters() { return _wavepacket_parameters; }
void MatlabData::wavepacket_parameters(WavepacketParameters *params) 
{ insist(params && !_wavepacket_parameters); _wavepacket_parameters = params; }

// CRP parameters
CRPParameters *MatlabData::crp_parameters() { return _crp_parameters; }
void MatlabData::crp_parameters(CRPParameters *params) 
{ insist(params && !_crp_parameters); _crp_parameters = params; }

// check Matlab data
void MatlabData::check_data()
{
  insist(MatlabData::r1()->n%2 == 0);
  insist(MatlabData::r2()->n%2 == 0);
  insist(MatlabData::theta()->n > MatlabData::wavepacket_parameters()->l_max);
  
  if(MatlabData::options()->calculate_reaction_probabilities) 
    insist(MatlabData::crp_parameters());
}

// destroy all data 

#define _FREE_(x) if(x) { delete x; x = 0; }

void MatlabData::destroy_all_data()
{
  std::cout << " Destroy Matlab data" << std::endl;
  
  _FREE_(_r1);
  _FREE_(_r2);
  _FREE_(_theta);
  _FREE_(_time);
  _FREE_(_options);
  _FREE_(_wavepacket_parameters);
  _FREE_(_crp_parameters);
  _FREE_(_si_coefficients);
  _potential = 0;
}

#ifdef _FREE_
#undef _FREE_
#endif
