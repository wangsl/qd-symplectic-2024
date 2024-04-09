
% $Id$

function [ V ] = H3PES(r1, r2, r3)

global H3PESName

if strcmp(H3PESName, 'LSTH')
  V = H3PESLSTH(r1, r2, r3);
elseif strcmp(H3PESName, 'LSTHFortran')
  V = H3PESLSTHFortran(r1, r2, r3);
elseif strcmp(H3PESName, 'BKMP2')
  V = H3PESBKMP2(r1, r2, r3);
else
  error('H3 potential energy surface name error');
end



