  
% $Id$

%function [] = main(jRot, nVib)

%close all
clear all
%clc

format long

%if nargin == 0 
  jRot = 0;
  nVib = 0;
%end

global H3PESName
global H2eV 
global H3Data

addpath('/home/wang/matlab/qd-cuda-symplectic-2018/H3');
addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/H3');
addpath('/home/wang/matlab/qd-cuda-symplectic-2018/common-cuda');
addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/common-cuda');
addpath('/home/wang/matlab/qd-cuda-symplectic-2018/common');
addpath('/home/wang/matlab/qd-cuda-symplectic-2018/build/common');

H3PESName = 'BKMP2';
%H3PESName = 'LSTH';
%H3PESName = 'LSTHFortran';

H2eV = 27.21138505;
MassAU = 1.822888484929367e+03;

dimensions = 2;

mH = 1.007825;
mD = 2.01410178;
mT = 3.0160492;

masses = [ mH, mH, mH ];

masses = masses*MassAU;

% time

time.total_steps = int32(20000000);
time.time_step = 0.1;
time.steps = int32(0);

% r1: R

r1.n = int32(512);
r1.r = linspace(0.4, 14.0, r1.n);
r1.left = r1.r(1);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = masses(1)*(masses(2)+masses(3))/(masses(1)+masses(2)+ ...
					   masses(3));

r1.dump = WoodsSaxon(4.0, 12.50, r1.r);

r1.r0 = 8.5;
r1.k0 = 12.0;
r1.delta = 0.12;

eGT = 1/(2*r1.mass)*(r1.k0^2 + 1/(2*r1.delta^2));
fprintf(' Gaussian wavepacket kinetic energy: %.15f\n', eGT)

% r2: r

r2.n = int32(512);
r2.r = linspace(0.4, 14.0, r2.n);
r2.left = r2.r(1);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = masses(2)*masses(3)/(masses(2)+masses(3));

r2.dump = WoodsSaxon(4.0, 12.5, r2.r);

% dividing surface

rd = 8.0;
nDivdSurf = int32((rd - min(r2.r))/r2.dr);
r2Div = double(nDivdSurf)*r2.dr + min(r2.r);
fprintf(' Dviding surface: %.8f\n', r2Div);

J = 0;
parity = 0;

% angle:
if dimensions == 2 
  theta.n = int32(1);
  theta.x = 1.0-1.0e-10;
  theta.w = 2.0;
  lMax = 0;
else
  theta.n = int32(120);
  [ theta.x, theta.w ] = GaussLegendreGridsMex(theta.n);
  lMax = 96;
end

options.wave_to_matlab = 'C2Matlab';
options.CRPMatFile = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
options.steps_to_copy_psi_from_device_to_host = int32(1000);
options.potential_cutoff = 0.2;
options.rotational_states = int32(0);
options.calculate_reaction_probabilities = int32(1);

potential = H3PESJacobi(r1.r, r2.r, theta.x, masses);

%J = 0;
%parity = 0;
%lMax = 120;

wavepacket_parameters.J = int32(J);
wavepacket_parameters.parity = int32(parity);
wavepacket_parameters.lMax = int32(lMax);

[ OmegaMin, OmegaMax ] = OmegaRange(J, parity, lMax);

wavepacket_parameters.OmegaMin = int32(OmegaMin);
wavepacket_parameters.OmegaMax = int32(OmegaMax);

[ psi, eH2 ] = InitWavePacket(r1, r2, theta, jRot, nVib);

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

% Reaction probabilities

CRP.mat_file = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
CRP.eDiatomic = eH2;
CRP.n_dividing_surface = nDivdSurf;
CRP.n_gradient_points = int32(51);
CRP.n_energies = int32(500);
eLeft = 0.3/H2eV;
eRight = 4.0/H2eV;
CRP.energies = linspace(eLeft, eRight, CRP.n_energies);
CRP.eta_sq = EtaSq(r1, CRP.energies-eH2);
CRP.CRP = zeros(size(CRP.energies));

% pack data to one structure

H3Data.r1 = r1;
H3Data.r2 = r2;
H3Data.theta = theta;
H3Data.potential = potential;
H3Data.time = time;
H3Data.options = options;
H3Data.CRP = CRP;

H3Data.wavepacket_parameters = wavepacket_parameters;

PlotPotWave();

% time evolution

tic
cudaSymplectic(H3Data);
toc



