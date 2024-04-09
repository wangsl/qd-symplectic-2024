
function [ p ] = fcnpak4b(nu1, nu2, mu, x)
assert(min(x) > -1, 'x should be > -1')
assert(max(x) < 1, 'x should be < 1')
p = fcnpakMex(int32(1), int32(nu1), int32(nu2), int32(mu), x);
p = p';

