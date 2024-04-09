
function [ psi, eH2, psiH2 ] = InitWavePacket(R1, R2, Theta, jRot, nVib)

n1 = R1.n;
n2 = R2.n;
nTheta = Theta.n;

r1 = R1.r;
delta = R1.delta;
r10 = R1.r0;
k0 = R1.k0;

G = (1/(pi*delta^2))^(1/4) * ...
    exp(-(r1-r10).^2/(2*delta*delta) - j*k0*r1);

fprintf(' Gaussian wavepacket module: %.15f\n', sum(conj(G).*G)*R1.dr);

[ eH2, phiH2 ] = H2VibRotWaveFunction(R2, jRot, nVib);

P = legendre(jRot, Theta.x, 'norm');

fprintf(' Legendre polynomail module: %.15f\n', ...
	sum(P.^2.*Theta.w));

GPhi = G'*phiH2;

GPhi = reshape(GPhi, [1, numel(GPhi)]);

psi = GPhi'*P;

psi = reshape(psi, [n1, n2, nTheta]);

