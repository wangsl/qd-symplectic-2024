
function [ V ] = DMBEIVPES(R1, R2, R3)

vOOMin = -0.19157004525;

V = DMBEIVMex(reshape(R1, [1, numel(R1)]), ...
	      reshape(R2, [1, numel(R2)]), ...
	      reshape(R3, [1, numel(R3)])) ...
    - vOOMin;

V = reshape(V, size(R1));


