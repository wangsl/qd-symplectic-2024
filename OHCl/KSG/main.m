

%function [] = main(jRot, nVib)

close all
clear all
clc

format long

%if nargin == 0 
jRot = 0;
nVib = 0;
%end

global H2eV 
global OHClData

fileID = fopen('/beegfs/work/wang/quantum-dynamics/m8n4.txt', 'r');
a = fscanf(fileID, '%f,', [8,1])

a = reshape(a, [1, length(a)]);
SI_coefficients.m = int32(length(a));
SI_coefficients.a = a;
SI_coefficients.b = fliplr(SI_coefficients.a);

H2eV = 27.21138505;

base_dir = '/home/wang/matlab/qd-cuda-symplectic-2018';

addpath(strcat(base_dir, '/common'));
addpath(strcat(base_dir, '/build/common'));
addpath(strcat(base_dir, '/common-cuda'));
addpath(strcat(base_dir, '/build/common-cuda'));
addpath(strcat(base_dir, '/OHCl'));
addpath(strcat(base_dir, '/build/OHCl'));

MassAU = 1.822888484929367e+03;

mH = 1.0079;
mCl = 35.453;
mO = 15.999;

masses = [ mO, mH, mCl ];

masses = masses*MassAU;

% time

time.total_steps = int32(20000000);
time.time_step = 1.0;
time.steps = int32(0);

% r1: R

r1.n = int32(768);
r1.r = linspace(0.4, 17.5, r1.n);
r1.left = r1.r(1);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = masses(1)*(masses(2)+masses(3))/(masses(1)+masses(2)+ masses(3));
r1.dump = WoodsSaxon(4.0, 16.0, r1.r);

r1.r0 = 9.5;
r1.k0 = 12.0;
r1.delta = 0.05;

eGT = 1/(2*r1.mass)*(r1.k0^2 + 1/(2*r1.delta^2));
fprintf(' Gaussian wavepacket kinetic energy: %.15f\n', eGT);

% r2: r

r2.n = int32(512);
r2.r = linspace(1.2, 14.0, r2.n);
r2.left = r2.r(1);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = masses(2)*masses(3)/(masses(2)+masses(3));
r2.dump = WoodsSaxon(4.0, 12.5, r2.r);

% dividing surface

rd = 9.0;
nDivdSurf = int32((rd - min(r2.r))/r2.dr);
r2Div = double(nDivdSurf)*r2.dr + min(r2.r);
fprintf(' Dviding surface: %.8f\n', r2Div);

% theta

theta.n = int32(180);
[ theta.x, theta.w ] = GaussLegendreGridsMex(theta.n);

% options

options.wave_to_matlab = 'OHClMatlab.m';
options.CRPMatFile = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
options.steps_to_copy_psi_from_device_to_host = int32(100);
options.potential_cutoff = 0.4;
options.rotational_states = int32(0);
options.calculate_reaction_probabilities = int32(1);

% setup potential energy surface and initial wavepacket
potential = OHClPESJacobi(r1.r, r2.r, theta.x, masses);

J = 0;
parity = 0;
lMax = 120;

wavepacket_parameters.J = int32(J);
wavepacket_parameters.parity = int32(parity);
wavepacket_parameters.lMax = int32(lMax);

[ OmegaMin, OmegaMax ] = OmegaRange(J, parity, lMax);

wavepacket_parameters.OmegaMin = int32(OmegaMin);
wavepacket_parameters.OmegaMax = int32(OmegaMax);

[ psi, eHCl ] = InitWavePacket(r1, r2, theta, jRot, OmegaMin, nVib);

P = AssociatedLegendreP(OmegaMin, OmegaMax, lMax, theta.x);
for k = 1 : theta.n
  P(k,:,:) = P(k,:,:)*sqrt(theta.w(k));
end

wavepacket_parameters.weighted_associated_legendres = P;

nOmegas = OmegaMax - OmegaMin + 1;
wavepackets = zeros([size(psi), nOmegas]);
wavepackets(:,:,:,1) = psi;

for k = 1 : theta.n
  wavepackets(:,:,k,:) = wavepackets(:,:,k,:)*sqrt(theta.w(k));
end

wavepacket_parameters.weighted_wavepackets = wavepackets;

% cummulative reaction probabilities

CRP.mat_file = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
CRP.eDiatomic = eHCl;
CRP.n_dividing_surface = nDivdSurf;
CRP.n_gradient_points = int32(51);
CRP.n_energies = int32(500);
eLeft = 0.01/H2eV + eHCl;
eRight = 1.0/H2eV + eHCl;
CRP.energies = linspace(eLeft, eRight, CRP.n_energies);
CRP.eta_sq = EtaSq(r1, CRP.energies-eHCl);
CRP.CRP = zeros(size(CRP.energies));

% pack data to one structure

OHClData.r1 = r1;
OHClData.r2 = r2;
OHClData.theta = theta;
OHClData.potential = potential;
OHClData.time = time;
OHClData.options = options;
OHClData.CRP = CRP;

HO2Data.SI_coefficients = SI_coefficients;

OHClData.wavepacket_parameters = wavepacket_parameters;

PlotPotWave();

% time evolution
tic
cudaSymplectic(OHClData);
toc

return

