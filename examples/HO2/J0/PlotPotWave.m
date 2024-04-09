
% $Id$

function [] = PlotPotWave()

global HO2Data

persistent has_PotWavePlot
persistent hpsi

persistent fig_index

theta = HO2Data.theta;

r1 = HO2Data.r1;
r2 = HO2Data.r2;
pot = HO2Data.potential;

omega = 1;

k = 1;

psiReal = HO2Data.wavepacket_parameters.weighted_wavepackets.real(:, :, k, omega)/sqrt(theta.w(k));
psiReal = psiReal';

if isempty(has_PotWavePlot)
  
  has_PotWavePlot = 1;
  
  figure(1);
  
  [ ~, hPES ] = contour(r1.r, r2.r, pot(:,:,k)', [ -0.2:0.01:0.3 ]);
  set(hPES, 'LineWidth', 0.75);
  set(hPES, 'LineColor', 'black');
  hold on;
  
  [ ~, hpsi ] = contour(r1.r, r2.r, psiReal, ...
			[ -3.0:0.02:-0.01 0.01:0.02:3.0 ], 'zDataSource', 'psiReal');
  set(hpsi, 'LineWidth', 1.0);
  set(gca, 'CLim', [-0.5, 0.5]);
  colormap jet
  colorbar('vert')
  hold off;
end

refreshdata(hpsi, 'caller');
drawnow

if isempty(fig_index) 
  fig_index = 0;
else
  fig_index = fig_index + 1;
end

index = strcat('000000000', int2str(fig_index));
fig_name = strcat('potwave/', index(end-8:end), '.png');

fprintf(' Figure: %s\n', fig_name)

figure(1);
print(fig_name, '-dpng');

