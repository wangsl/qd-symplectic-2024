
% J. Chem. Phys. 108, 5404 (1998)
% for J = 0

function [] = HO2main(jRot, nVib)

%clear all
%close all
%clc

format long

[ ~, host ] = system('hostname');
fprintf(' Host: %s\n', strtrim(host));
fprintf(' Current Matlab PID: %d\n', feature('getpid'));
fprintf(' Working directory: %s\n', pwd);
fprintf('\n');

baseFolder = '/scratch/work/wang/centos-8/matlab/qd-cuda-symplectic-2018/';

addpath([baseFolder, 'common']);
addpath([baseFolder, 'build/common']);
addpath([baseFolder, 'common-cuda']);
addpath([baseFolder, 'build/common-cuda']);

addpath([baseFolder, 'HO2']);
addpath([baseFolder, 'build/HO2']);

if nargin == 0 
  jRot = 1;
  nVib = 0;
end

global H2eV 
global HO2Data

H2eV = 27.21138505;

MassAU = 1.822888484929367e+03;

mH = 1.0079;
mO = 15.999;
 
masses = [ mH, mO, mO ];

masses = masses*MassAU;

% time

time.total_steps = int32(2000000);
time.time_step = 0.5;
time.steps = int32(0);

% r1: R

r1.n = int32(512);
r1.r = linspace(0.6, 16.0, r1.n);
r1.left = r1.r(1);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = masses(1)*(masses(2)+masses(3))/(masses(1)+masses(2)+ masses(3));
r1.dump = WoodsSaxon(4.0, 14.50, r1.r);

r1.r0 = 8.5;
r1.k0 = 13;
r1.delta = 0.6;

eGT = 1/(2*r1.mass)*(r1.k0^2 + 1/(2*r1.delta^2));
fprintf(' Gaussian wavepacket kinetic energy: %.15f\n', eGT)

% r2: r

r2.n = int32(512);
r2.r = linspace(1.3, 18.0, r2.n);
r2.left = r2.r(1);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = masses(2)*masses(3)/(masses(2)+masses(3));

r2.dump = WoodsSaxon(4.0, 16.5, r2.r);

% dividing surface

rd = 10.4;
nDivdSurf = int32((rd - min(r2.r))/r2.dr);
r2Div = double(nDivdSurf)*r2.dr + min(r2.r);
fprintf(' Dviding surface: %.15f\n', r2Div);

% theta

theta.n = int32(180);
[ theta.x, theta.w ] = GaussLegendreGridsMex(theta.n);

% options

options.wave_to_matlab = 'HO2Matlab';
options.CRPMatFile = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
options.steps_to_copy_psi_from_device_to_host = int32(100);
options.potential_cutoff = 0.5;
options.rotational_states = int32(1);
options.calculate_reaction_probabilities = int32(1);
options.catastrophe_criterion_with_module = 1.1;
options.converged_criterion_with_module = 1e-4;

% setup potential energy surface and initial wavepacket
potential = DMBEIVPESJacobi(r1.r, r2.r, theta.x, masses);

J = 0;
parity = 0;
lMax = 120;

wavepacket_parameters.J = int32(J);
wavepacket_parameters.parity = int32(parity);
wavepacket_parameters.lMax = int32(lMax);

[ OmegaMin, OmegaMax ] = OmegaRange(J, parity, lMax);

wavepacket_parameters.OmegaMin = int32(OmegaMin);
wavepacket_parameters.OmegaMax = int32(OmegaMax);

[ psi, eO2 ] = InitWavePacket(r1, r2, theta, jRot, OmegaMin, nVib);

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

wavepacket_parameters.weighted_wavepackets.real = real(wavepackets);
wavepacket_parameters.weighted_wavepackets.imag = imag(wavepackets);

% Reaction probabilities

CRP.mat_file = sprintf('CRPMat-j%d-v%d.mat', jRot, nVib);
CRP.eDiatomic = eO2;
CRP.n_dividing_surface = nDivdSurf;
CRP.n_gradient_points = int32(51);
CRP.n_energies = int32(500);
eLeft = 0.8/H2eV;
eRight = 1.8/H2eV;
CRP.energies = linspace(eLeft, eRight, CRP.n_energies);
CRP.eta_sq = EtaSq(r1, CRP.energies-eO2);
CRP.CRP = zeros(size(CRP.energies));

% pack data to one structure

HO2Data.r1 = r1;
HO2Data.r2 = r2;
HO2Data.theta = theta;
HO2Data.potential = potential;
HO2Data.time = time;
HO2Data.options = options;
HO2Data.CRP = CRP;

HO2Data.wavepacket_parameters = wavepacket_parameters;

PlotPotWave();

% time evolution

tic
cudaSymplectic(HO2Data);
fprintf(' ');
toc

fprintf(' Simulation has finished\n');

return

