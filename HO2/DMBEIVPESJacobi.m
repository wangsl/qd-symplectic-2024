
% $Id$

function [ V ] = DMBEIVPESJacobi(R, r, cosTheta, masses)

f2 = masses(3)/(masses(2)+masses(3));
f3 = 1 - f2;

[ y, x, z ] = meshgrid(r, R, cosTheta);

r23 = y;
r13 = sqrt(x.^2 + (f3*y).^2 - 2*f3*x.*y.*z);
r12 = sqrt(x.^2 + (f2*y).^2 + 2*f2*x.*y.*z);

V = DMBEIVPES(r23, r12, r13);

