
function [ V ] = OOPES(rOO)

% rOO = R1
% angle: OOH=0.0 R2=R1+R3
% angle: OOH=180.0 R3=R1+R2

%r1 = reshape(rOO, [1, numel(rOO)]);

r1 = rOO;

%r3 = zeros(size(r1));
%r3(:) = 100.0;

r3 = 100.0*ones(size(r1));

r2 = r1 + r3;

vOOMin = -0.19157004525;

V = DMBEIVMex(r1, r2, r3) - vOOMin;
