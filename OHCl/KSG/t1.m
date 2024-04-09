
format long

addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/OHCl');

HClPES(1.2835)

r = linspace(1.0, 10.0, 128);

v = HClPES(r);

plot(r, v)

