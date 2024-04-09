

global H3PESName

a0 = 1/0.5291772049;

%r1 = 100;
%r2 = r;
%r3 = r1+r2;

H3PESName = 'LSTHFortran';
[ r, v ] = fminsearch(@(r)H3PES(100, r, 100+r), 2.5)

H3PESName = 'LSTH';
[ r, v ] = fminsearch(@(r)H3PES(100, r, 100+r), 2.5)
%H3PES(r1, r2, r3)

H3PESName = 'BKMP2';
[ r, v ] = fminsearch(@(r)H3PES(100, r, 100+r), 2.5)
%H3PES(r1, r2, r3)

