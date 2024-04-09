
function [ V ] = OHClRKHSPES(rOH, rOCl, rHCl)

V = OHClRKHSMex(reshape(rOH, [1, numel(rOH)]), ...
		reshape(rOCl, [1, numel(rOCl)]), ...
		reshape(rHCl, [1, numel(rHCl)]));

V = reshape(V, size(rOH));


