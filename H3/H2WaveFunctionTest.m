
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

addpath(genpath('/home/wang/matlab/quantum-dynamics'));

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

MassAU = 1.822888484929367e+03;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

r.n = int32(2048);
r.r = linspace(0.4, 20.0, r.n);
r.dr = r.r(2) - r.r(1);
r.mass = mH/2;

energy_levels = 0:1:15; %[ 0 1 6 13 14 15 16 ];

[ e, psi ] = H2VibRotWaveFunction(r, 0, energy_levels);

plot(r.r, psi, 'LineWidth', 2);
legend(arrayfun(@num2str, energy_levels, 'UniformOutput', false));

print -dpdf H2-wavefunction.pdf

e + vH2Min

return

  -0.164566592182935
  -0.145607293917483
  -0.127720805943968
  -0.110884360261777
  -0.095083828470802
  -0.080314527140222
  -0.066582402807662
  -0.053905749601210
  -0.042317668977451
  -0.031869560175718
  -0.022636068473634
  -0.014722194285049
  -0.008273764988278
  -0.003493015452886
  -0.000659266463932
   0.000036144518170