
% $Id$

function [ e, phi ] = OOVibRotWaveFunction(R, jRot, varargin)

%% n should be converted to double, otherwise it will be wrong

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

V = OOPES(r) + jRot*(jRot+1)./(2*mu*r.^2);

I = 1:n;
V = V + 1/(2*mu*dr*dr)*(pi^2/3 - 1./(2*I.^2));

% set diagonal elements
H(1:n+1:end) = V;

[vecs, energies] = eig(H);

nVbs = nVbs+1;
e = diag(energies);
e = e(nVbs);
phi = vecs(:, nVbs)/sqrt(dr);

phi = phi';

ePot = sum(phi.*OOPES(r).*phi)*dr;
eRot = sum(phi.*jRot*(jRot+1)./(2*mu*r.^2).*phi)*dr;
eKin = e - ePot - eRot;

fprintf(' O2 vibrational wavefunction module: %.15f\n', ...
	sum(phi.*phi)*dr);
fprintf(' O2 total energy: %.15f\n', e);
fprintf(' O2 radial kinetic energy: %.15f\n', eKin);
fprintf(' O2 rotational energy: %.15f\n', eRot);
fprintf(' O2 potential energy: %.15f\n', ePot);

