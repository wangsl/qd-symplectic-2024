
% $Id$

function [ e, phi ] = HClVibRotWaveFunction(R, jRot, varargin)

n = double(R.n);
dr = R.dr;
mu = R.mass;
r = R.r;

if length(varargin) > 0 
  nVbs = [ varargin{:} ];
else
  nVbs = 0:1:n-1;
end

[I, J] = ind2sub([n, n], 1:n^2);

H = 1/(2*mu*dr*dr)*(-1).^(I-J).*(2./(I-J).^2 - 2./(I+J).^2);

H = reshape(H, [n, n]);

V = HClPES(r) + jRot*(jRot+1)./(2*mu*r.^2);

I = 1:n;
V = V + 1/(2*mu*dr*dr)*(pi^2/3 - 1./(2*I.^2));

% set diagonal elements
H(1:n+1:end) = V;

[ vecs, energies ] = eig(H);

nVbs = nVbs+1;
e = diag(energies);
e = e(nVbs);
phi = vecs(:, nVbs)/sqrt(dr);

phi = phi';

ePot = sum(phi.*HClPES(r).*phi)*dr;
eRot = sum(phi.*jRot*(jRot+1)./(2*mu*r.^2).*phi)*dr;
eKin = e - ePot - eRot;

fprintf(' HCl vibrational wavefunction module: %.15f\n', ...
        sum(phi.*phi)*dr);
fprintf(' HCl total energy: %.15f\n', e);
fprintf(' HCl radial kinetic energy: %.15f\n', eKin);
fprintf(' HCl rotational energy: %.15f\n', eRot);
fprintf(' HCl potential energy: %.15f\n', ePot);



