
function [ V ] = Cl2PES(r)

R2 = r;
R1 = zeros(size(R2));
R1(:) = 100.0;

V = HCl2GHNS(R1, R2, R1+R2);

VCl2Min = -0.09067089035558364815; %2.4673 eV;
V = V - VCl2Min;

return
