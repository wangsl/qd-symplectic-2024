
% $Id$

% Ref: J. Chem. Phys. 104. 7139 (1996)

function [ varargout ] = H3PESLSTHFortran(r1, r2, r3)

persistent first

if isempty(first) 
  fprintf(' To use LSTH Fortran PES\n');
  first = 0;
end

vH2Min = -1.1744744 + 1;

varargout{1} = H3LSTHMex(r1, r2, r3);

varargout{1} = varargout{1} - vH2Min;


