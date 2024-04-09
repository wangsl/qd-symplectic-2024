
clear all
close all
clc

format long

addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/OHCl');

MassAU = 1.822888484929367e+03;

mH = 1.0079;
mCl = 35.453;
mO = 15.999;

masses = [ mO, mCl, mH ];

masses = masses*MassAU;

r2.n = int32(16384);
r2.r = linspace(1.0, 40.0, r2.n);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = masses(2)*masses(3)/(masses(2)+masses(3));

energy_levels = 0:1:28; %[ 0 1 6 13 14 15 16 ];
energy_levels = [ 0 ];

[ e, psi ] = HClVibRotWaveFunction(r2, 0, energy_levels);

plot(r2.r, psi, 'LineWidth', 2);
legend(arrayfun(@num2str, energy_levels, 'UniformOutput', false));


return

  -0.162955955345964
  -0.149878080963364
  -0.137347159814628
  -0.125363191899771
  -0.113926177218793
  -0.103036115771709
  -0.092693007558507
  -0.082896852579162
  -0.073647650833726
  -0.064945402322157
  -0.056790107044470
  -0.049181765000677
  -0.042120376190745
  -0.035605940614699
  -0.029638458272542
  -0.024217929164249
  -0.019344353289864
  -0.015017730649347
  -0.011238061242712
  -0.008005345069973
  -0.005319582131093
  -0.003180772426095
  -0.001588915954974
  -0.000544012717731
  -0.000046062714407
