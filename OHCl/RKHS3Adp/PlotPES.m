
close all;
clear all;
clc

addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/OHCl/RKHS3Adp')

mH = 1.0079;
mCl = 35.453;
mO = 15.999;

MassAU = 1.822888484929367e+03;

masses = [ mO, mCl, mH ];
masses = masses*MassAU;

r = linspace(1.4, 9, 128);
R = linspace(1.2, 12, 128);

set(gcf, 'PaperOrientation', 'landscape');

n = 3;
m = 3;
k = 0;

gamma0 = 0.001/180*pi;

for i = 1 : n
  for j = 1 : m
    k = k + 1
    subplot(n, m, k)
    gamma = 20.0/180*pi*(k-1)+gamma0
    V = OHClPESJacobi(R, r, cos(gamma), masses);
    [C, h] = contour(R, r, V',  [-0.2:0.01:0.4]);
    set(h, 'LineWidth', 0.75);
    set(h, 'LineColor', 'black');
    axis normal;
  end
end
    
print -dpdf PES-gamma.pdf
