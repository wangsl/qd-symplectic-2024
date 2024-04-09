
%close all;
clear all;
clc

addpath('/home/wang/matlab/quantum-dynamics/build/OHCl')

mH = 1.0079;
mCl = 35.453;
mO = 15.999;

masses = [ mO, mH, mCl ];

r = linspace(1.8, 9, 256);
R = linspace(1.8, 12, 256);

set(gcf, 'PaperOrientation', 'landscape');

n = 3;
m = 3;
k = 0;

for i = 1 : n
  for j = 1 : m
    k = k + 1
    subplot(n, m, k)
    gamma0 = 0.0;
    if k == 1 
      gamma0 = 0.1/180*pi;
    end
    V = OHClPESJacobi(R, r, 5.0/180*pi*(k-1)+gamma0, masses);
    [C, h] = contour(r, R, V,  [-0.2:0.01:0.2]);
    set(h, 'LineWidth', 0.75);
    set(h, 'LineColor', 'black');
    axis normal;
  end
end
    
print -dpdf PES-gamma.pdf
