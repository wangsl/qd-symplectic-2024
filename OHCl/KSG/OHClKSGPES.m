
function [ V ] = OHClKSGPES(rOH, rOCl, rHCl)

vHClMin = -0.1697;
vHClMin = 0.0;

V = OHClKSGMex(reshape(rOH, [1, numel(rOH)]), ...
	       reshape(rOCl, [1, numel(rOCl)]), ...
	       reshape(rHCl, [1, numel(rHCl)])) ...
    - vHClMin;

V = reshape(V, size(rOH));

