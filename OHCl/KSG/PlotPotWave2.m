
% $Id$

function [] = PlotPotWave2(r1, r2, pot)

%persistent has_PotWavePlot
%persistent hpsi

%[ psiReal, psiImag ] = WavePacket(psi);

k = 1;

figure(1);

[ ~, hPES ] = contour(r1.r, r2.r, pot(:,:,k)', [ -0.4:0.01:0.4 ]);
set(hPES, 'LineWidth', 0.75);
set(hPES, 'LineColor', 'black');
%axis square
hold on;

%refreshdata(hpsi, 'caller');
%drawnow

