
% $Id$

function [ V ] = H2PES(r)

R = 100.0*ones(size(r));

V = H3PES(r, R, r+R);

return
