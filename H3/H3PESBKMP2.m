
% $Id$

% Ref: J. Chem. Phys. 104. 7139 (1996)

function [ varargout ] = H3PESBKMP2(r1, r2, r3)

persistent first

if isempty(first) 
  fprintf(' To use BKMP2 PES\n');
  first = 0;
end

vH2Min = -0.174495770896975;

if nargout == 0 | nargout == 1
  varargout{1} = BKMP2Mex(r1, r2, r3);
  %elseif nargout == 2
  %  [ varargout{1}, varargout{2} ] = BKMP2Mex(R1, R2, R3);
end

varargout{1} = varargout{1} - vH2Min;
