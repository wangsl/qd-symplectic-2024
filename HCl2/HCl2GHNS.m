
%===============================================
% the analytical fit form of the HCl2 potential
% Miguel Gonzalez, J.Hijazo, J.J.Novoa,and R.Sayos
% J. Chem. Phys 108, 3168 (1998)
%===============================================
%
% there are some mistakes in the original paper,
% instead of (rho(1),rho(2),rho(3)),
% the following internal coordinates are used
% in eq(5), both polynomials eq(6)
%           and range function eq(7)
% s1=0.21007*rho(1)+0.95485*rho(2)+0.21007*rho(3)
% s2=0.70711*rho(1)-0.70711*rho(3)
% s3=0.67518*rho(1)+0.29708*rho(2)-0.67518*rho(3)
% those information is from private communication
% E-mail: miguel@dymamics.qf.ub.es
% (from the bottom of the first page in the paper)

% R1 = RHCl
% R2 = RClCl
% R3 = RHCl

function [ VHCl2 ] = HCl2GHNS(R1, R2, R3)

% convert Rs from Bohr to Angstrom
B2A = 0.529177249;
R1 = R1*B2A;
R2 = R2*B2A;
R3 = R3*B2A;

VHCl2 = VHCl(R1) + VCl2(R2) + VHCl(R3) + V3HCl2(R1, R2, R3);

% convert potential energy from eV to atomic units
H2eV = 27.21138505;
VHCl2 = VHCl2/H2eV;
return

function [ V ] = VCl2(R)
De = 2.4673;
Re = 1.9939;
a1 = 5.2714;
a2 = 9.5082;
a3 = 9.5665;
a4 = 13.0256;
a5 = 18.4929;
rho = R - Re;
V = -De*(1+rho.*(a1+rho.*(a2+rho.*(a3+rho.*(a4+rho.*a5))))).*exp(-a1*rho);
return


function [ V ] = VHCl(R)
De = 4.6123;
Re = 1.2687;
a1 = 3.6314;
a2 = 2.9360;
a3 = 1.0346;
a4 = -0.8237;
a5 = 1.6907;
rho = R- Re;
V = -De*(1+rho.*(a1+rho.*(a2+rho.*(a3+rho.*(a4+rho.*a5))))).*exp(-a1*rho);
return


function [ P ] = PHCl2(s1, s2, s3)

V0 = 0.5840;

c = [  0.42710 -0.11400  0.82660 -0.54600  0.11980 ...
       0.00000  0.00000  0.00000  0.00000  0.00000 ...
       0.89540  0.69570  0.50000 -0.68770  0.00000 ...
       0.00000  0.00000 -0.02378 -0.39740  0.00000 ...
       2.50870 -2.71150 -0.55000  1.19860 -0.14970 ...
       0.00000  0.00000  0.00000  0.00000  0.83020 ...
       0.44740 -0.47650  0.00000  0.00000 -0.07957 ...
       2.15540 -16.3816 13.07980 -2.10400  0.00000 ...
       0.00000  0.00000  3.01820 -0.98590  0.00000 ...
      -8.20580 -0.07202  2.74710  0.00000  0.00000 ...
       1.28580 -4.47950  3.58280  0.00000 -0.40510 ];

s = 1.0;
ii = 0;
for i = 0 : 5 
  s1i = s1.^i;
  for j = 0 : 5-i 
    s2j = s2.^j;
    for k = 0 : 5-i-j
      if(i+j+k > 0)
	ii = ii+1;
	s3k = s3.^k;
	s = s + c(ii).*s1i.*s2j.*s3k;
      end
    end
  end
end

P = V0*s;
return


function [ T ] = THCl2(s1, s2, s3)
gamma1 = 3.4515;
gamma2 = 0.0;
gamma3 = -2.2231;
T = (1-tanh(gamma1*s1/2)) .* ...
    (1-tanh(gamma2*s2/2)) .* ...
    (1-tanh(gamma3*s3/2));
return


function [ V ] = V3HCl2(R1, R2, R3)
R10 = 3.1157;
R20 = 2.0035;
R30 = 3.1157;
rho1 = R1 - R10;
rho2 = R2 - R20;
rho3 = R3 - R30;
s1 = 0.21007*rho1 + 0.95485*rho2 + 0.21007*rho3;
s2 = 0.70711*rho1 - 0.70711*rho3;
s3 = -0.67518*rho1 + 0.29708*rho2 - 0.67518*rho3;
V = PHCl2(s1, s2, s3) .* THCl2(s1, s2, s3);
return

