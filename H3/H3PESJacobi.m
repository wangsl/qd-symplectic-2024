
% $Id$

function [ V ] = H3PESJacobi(R, r, cosTheta, masses)

global H3PESName

f2 = masses(3)/(masses(2)+masses(3));
f3 = 1 - f2;

if strcmp(H3PESName, 'LSTH') | strcmp(H3PESName, 'LSTHFortran')
  rMin = 0.4;
  VMax = 2.5;
elseif strcmp(H3PESName, 'BKMP2')
  rMin = 0.2;
  VMax = 2.5;
else
  error('H3 potential energy surface name error');
end

nR = length(R);
nr = length(r);
nTheta = length(cosTheta);

V = zeros(nR, nr, nTheta);
V(:) = VMax;

[ y, x, z ] = meshgrid(r, R, cosTheta);

r23 = y;
r13 = sqrt(x.^2 + (f3*y).^2 - 2*f3*x.*y.*z);
r12 = sqrt(x.^2 + (f2*y).^2 + 2*f2*x.*y.*z);

grids = find(r12>rMin & r13>rMin & r23>rMin);

V(grids) = H3PES(r12(grids), r23(grids), r13(grids));



