
function [ V ] = HClPES(rHCl)

% can not be linear

rOH = 100*ones(size(rHCl));
rOCl = rHCl + rOH - 0.0001;

V = OHClKSGPES(rOH, rOCl, rHCl);

