
/* created at: 2021-05-18 19:57:41 */

#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "matlabStructures.h"
#include "die.h"

ostream & operator <<(ostream &s, const RadialCoordinate &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void RadialCoordinate::write_fields(ostream &s) const
{
  s << Indent() << "n " << n << "\n";
  s << Indent() << "left " << left << "\n";
  s << Indent() << "dr " << dr << "\n";
  s << Indent() << "mass " << mass << "\n";
}

ostream & operator <<(ostream &s, const AngleCoordinate &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void AngleCoordinate::write_fields(ostream &s) const
{
  s << Indent() << "n " << n << "\n";
  s << Indent() << "x " << x << "\n";
  s << Indent() << "w " << w << "\n";
}

ostream & operator <<(ostream &s, const EvolutionTime &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void EvolutionTime::write_fields(ostream &s) const
{
  s << Indent() << "total_steps " << total_steps << "\n";
  s << Indent() << "time_step " << time_step << "\n";
  s << Indent() << "steps " << steps << "\n";
}

ostream & operator <<(ostream &s, const SICoefficients &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void SICoefficients::write_fields(ostream &s) const
{
  s << Indent() << "a " << a << "\n";
  s << Indent() << "b " << b << "\n";
}

ostream & operator <<(ostream &s, const Options &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void Options::write_fields(ostream &s) const
{
  if (wave_to_matlab)
    s << Indent() << "wave_to_matlab " << wave_to_matlab << "\n";
  s << Indent() << "steps_to_copy_psi_from_device_to_host " << steps_to_copy_psi_from_device_to_host << "\n";
  s << Indent() << "potential_cutoff " << potential_cutoff << "\n";
  s << Indent() << "calculate_reaction_probabilities " << calculate_reaction_probabilities << "\n";
  s << Indent() << "rotational_states " << rotational_states << "\n";
  s << Indent() << "converged_criterion_with_module " << converged_criterion_with_module << "\n";
  s << Indent() << "catastrophe_criterion_with_module " << catastrophe_criterion_with_module << "\n";
}

ostream & operator <<(ostream &s, const WavepacketParameters &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void WavepacketParameters::write_fields(ostream &s) const
{
  s << Indent() << "J " << J << "\n";
  s << Indent() << "parity " << parity << "\n";
  s << Indent() << "l_max " << l_max << "\n";
  s << Indent() << "omega_min " << omega_min << "\n";
  s << Indent() << "omega_max " << omega_max << "\n";
}

ostream & operator <<(ostream &s, const CRPParameters &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void CRPParameters::write_fields(ostream &s) const
{
  s << Indent() << "n_dividing_surface " << n_dividing_surface << "\n";
  s << Indent() << "n_gradient_points " << n_gradient_points << "\n";
  s << Indent() << "n_energies " << n_energies << "\n";
}

