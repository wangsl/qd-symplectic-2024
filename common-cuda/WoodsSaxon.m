
% $Id$

function [ f ] = WoodsSaxon(Cd, xd, x)
f = 1./(1+exp(Cd*(x-xd)));
return

