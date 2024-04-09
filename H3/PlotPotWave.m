
% $Id$

function [] = PlotPotWave()

global H3Data

persistent has_PotWavePlot
persistent hpsi

theta = H3Data.theta;

r1 = H3Data.r1;
r2 = H3Data.r2;
pot = H3Data.potential;

omega = 1;
k = 1;

psiReal = real(H3Data.wavepacket_parameters.weighted_wavepackets(:, ...
						  :, k, omega))/sqrt(theta.w(k));
psiReal = psiReal';

if isempty(has_PotWavePlot)
  
  has_PotWavePlot = 1;
  
  figure(1);
  
  [ ~, hPES ] = contour(r1.r, r2.r, pot(:,:,k)', [ -0.2:0.01:0.3 ]);
  set(hPES, 'LineWidth', 0.75);
  set(hPES, 'LineColor', 'black');
  %axis square
  hold on;
  
  [ ~, hpsi ] = contour(r1.r, r2.r, psiReal, ...
			[ -2.0:0.02:-0.02 0.02:0.02:2.0 ], 'zDataSource', 'psiReal');
  set(hpsi, 'LineWidth', 1.5);
  set(gca, 'CLim', [-0.5, 0.5]);
  %axis square
  colormap jet
  colorbar vert
  hold off;
end

refreshdata(hpsi, 'caller');
drawnow






