
function [ OmegaMin, OmegaMax ] = OmegaRange(J, p, lMax)

OmegaMax = min(J, lMax);

if J == 0 
  assert(p == 0);
  OmegaMin = 0;
  OmegaMax = 0;
elseif rem(J+p, 2) == 0 
  OmegaMin = 0;
else
  OmegaMin = 1;
end

