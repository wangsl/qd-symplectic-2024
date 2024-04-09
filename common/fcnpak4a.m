
function [ p ] = fcnpak4a(nu, mu1, mu2, x)
assert(min(x) > -1, 'x should be > -1')
assert(max(x) < 1, 'x should be < 1')
p = fcnpakMex(int32(0), int32(nu), int32(mu1), int32(mu2), x);
p = p';

