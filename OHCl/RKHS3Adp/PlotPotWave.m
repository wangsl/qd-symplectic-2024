
% $Id$

function [] = PlotPotWave()

global OHClData

persistent has_PotWavePlot
persistent hpsi
persistent fig_index

theta = OHClData.theta;

r1 = OHClData.r1;
r2 = OHClData.r2;

omega = 1;

k = theta.n;

pot = OHClData.potential(:,:,k)';

psiReal = real(OHClData.wavepacket_parameters.weighted_wavepackets(:, ...
                                                  :, k, omega))/sqrt(theta.w(k));
psiReal = psiReal';

if isempty(has_PotWavePlot)
  
  has_PotWavePlot = 1;
  
  figure(1);
  
  [ ~, hPES ] = contour(r1.r, r2.r, pot, [ -0.1:0.01:0.4 ]);
  set(hPES, 'LineWidth', 0.75);
  set(hPES, 'LineColor', 'black');
  %axis square
  hold on;
  
  [ ~, hpsi ] = contour(r1.r, r2.r, psiReal, ...
			[ -2.0:0.02:-0.01 0.01:0.02:1.0 ], 'zDataSource', 'psiReal');
  set(hpsi, 'LineWidth', 1.5);
  set(gca, 'CLim', [-0.5, 0.5]);
  %axis square
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

index = strcat('00000000', int2str(fig_index));
fig_name = strcat('potwave/', index(end-8:end), '.png');

fprintf(' Figure: %s\n', fig_name)

figure(1);
print(fig_name, '-dpng');

