

function [ P ] = AssociatedLegendreP(OmegaMin, OmegaMax, lMax, x)

assert(OmegaMin == 0 | OmegaMin == 1);
assert(OmegaMax >= OmegaMin);
assert(lMax >= OmegaMax);
assert(length(x) >= lMax);

lMin = OmegaMin;

fprintf('\n Use fcnpak for associated Legendre polynomials\n')

P = zeros(numel(x), lMax-lMin+1, OmegaMax-OmegaMin+1);
for Omega = OmegaMin : OmegaMax
  P(:, :, Omega-OmegaMin+1) = fcnpak4b(lMin, lMax, Omega, x);
end

return

for l = lMin : lMax
  p = legendre(l, x, 'norm');
  OmegaStart = 1 + OmegaMin;
  OmegaEnd = min(OmegaMax+1, l+1);
  P(:, l-lMin+1, 1:OmegaEnd-OmegaStart+1) = p(OmegaStart:OmegaEnd, :)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

close all
clear all

n = 220;
 
[ x, w ] = GaussLegendreGrids(n);

OmegaMin = 1;
OmegaMax = 32;
lMax = 180;

tic
P = AssociatedLegendreP(OmegaMin, OmegaMax, lMax, x);
toc

for k = 1 : n
  P(k, :, :) = P(k, :, :)*sqrt(w(k));
end

Om = 32;

iOffset = 1 - OmegaMin;

Om = Om + iOffset;

for l1 = Om-4 : Om+5
  for l2 = Om-4 : Om+5
    fprintf(' %4d %4d  %16.12f\n', l1-iOffset, l2-iOffset, sum(P(:, l1, Om).*P(:, l2, Om)))
  end
end

%}