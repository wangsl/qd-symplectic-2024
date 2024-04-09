
% $Id$

function [ V ] = HCl2PESJacobi(R, r, Theta, masses)

f2 = masses(3)/(masses(2)+masses(3));
f3 = 1 - f2;

nR = length(R);
nr = length(r);
nTheta = length(Theta);

[x, y, z] = meshgrid(R, r, Theta);

r23 = y;
r12 = sqrt((f2*y).^2 + x.^2 - 2*f2*y.*x.*cos(z));
r13 = sqrt((f3*y).^2 + x.^2 + 2*f3*y.*x.*cos(z));

r23 = reshape(r23, [numel(r23), 1]);
r12 = reshape(r12, [numel(r12), 1]);
r13 = reshape(r13, [numel(r13), 1]);

% r12 : rHCl : R1
% r13 : rHCl : R3
% r23 : rClCl : R2

V = HCl2GHNS(r12, r23, r13);

V = reshape(V, [nr, nR, nTheta]);

V = permute(V, [2 1 3]);

VCl2Min = -0.09067089035558364815; %2.4673 eV;
V = V - VCl2Min;

return

