	subroutine prepot

c   RP SURFACES
c ----------------------------------------------------------------------------
c..						Version Date: December 17,2002
C ---DERIVATIVES NOT IMPLEMENTED

c.. The O(3P) + HCl --> OH + Cl MRCI+Q/CBS(aug-cc-pVnZ; n = 2-4) 3A'
c.. potential surface of B. Ramachandran and K. A. Peterson.

c.. The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
c.. to interpolate the surface consisting of 490 geometries spanning
c.. O-H-Cl angles of 70-180 deg.  This surface does not contain the
c.. H + ClO arrangement and no claims of accuracy is made for O-H-Cl
c.. angles smaller than 70 deg. 
c ----------------------------------------------------------------------------
c  USAGE:

c  All arguments are passed through the common block /potcm/. On input,
c  On Input:
c	r(1) = R(O-H)
c	r(2) = R(H-Cl)
c	r(3) = R(Cl-O).
c
c On Output:
c	energy* = Energy in Hartree w.r to the asymptotic O + HCl minimum.
C ---DERIVATIVES NOT YET IMPLEMENTED
c	dedr   = Derivatives of the potential w.r. to r(1), r(2), and r(3).

c  NOTE:
c    Before any actual potential energy calculations are made, a single
c    call to prepot must be made:
c	call prepot

c    This reads in the two body and three body expansion coefficients,
c..  computes the asymptotic O + HCl minimum energy and saves it. Later,
c    the potential energy is computed by calling pot:
c	call pot

c ----------------------------------------------------------------------------
	implicit double precision (a-h, o-z)
	parameter (n3=490)
	character*40 file1
	common /potcm/r(3), energy, dedr(3)
	common/geoms/R3(n3,3)
	common/kernel/nk, mk, ns
	common/cffts/C3(n3), e0, sf
	common/v2sum/R1g, R2g, R3g, VOH, VHCl, VClO
	data R1a,R2a,theta/50.0d0,2.41041209d0,180.0d0/
	data zero,one,two,eps/0.0d0,1.0d0,2.0d0,1.0d-12/

	file1 = 'threebody_3Ap_490.txt'

	open (2,file=file1,form='formatted',status='old')
	read (2,*) sf, nk, mk, ns
	do i = 1,n3
	    read (2,'(3f10.2,d25.15)') (R3(i,j),j=1,3), C3(i)	
	enddo

c.. Calculate the asymptotic O + HCl minimum energy e0:
c.. Sum of the two body terms:
	R1g = R1a
	R2g = R2a
	R3g = R1a+R2a
	call twobody

	e0 = VOH + VHCl + VClO

	write (6,*)
	write (6,*)
	write(6,*) 'PREPOT has been called for the OHCl RP 3A'' surface.'
	write(6,*) 'Version date December 17, 2002'
	write(6,*) 'Three-body coefficients:'
	write(6,'(a)') file1
	write(6,*)
	write (6,101) nk, mk, ns
	write (6,*) 'Asymptotic Potential: ',e0
	write (6,*)
	return

	entry pot

c.. Twobody terms:
	R1g = r(1)
	R2g = r(2)
	R3g = r(3)
	call twobody	

c.. Threebody term:
	xg = r(1)**ns
	yg = r(2)**ns
	zg = (one - CFG(r(1),r(2),r(3)))/two + eps
	energy = zero
	do j = 1,n3
            xj = R3(j,1)**ns
            yj = R3(j,2)**ns
	    zj = (one - cosd(R3(j,3)))/two + eps
            energy = energy + C3(j)*q(xg,xj)*q(yg,yj)*qK(zg,zj)
	enddo

	energy = sf*energy + VOH + VHCl + VClO - e0

	return
101	format (/10x,'Order of the kernel n = ',i3/10x,
     .	'Reciprocal power m = ',i3/10x,
     .	'Power of the r coordinate (r**s)  s = ',i3/)
	end
c ----------------------------------------------------------------------------
	subroutine twobody
	implicit double precision (a-h, o-z)
	common/v2sum/R1g, R2g, R3g, VOH, VHCl, VClO
	dimension C1(6), C2(6)

c.. SEC-scaled MRCI+Q/CBS OH potential:
	data (C1(i),i=1,6)/0.17046982d0,2.45278241d0,1.52360833d0,
     .	0.62204843d0,1.83107461d0,0.0000005274d0/

c.. Unscaled MRCI+Q/CBS HCl potential:
	data (C2(i),i=1,6)/0.17035320d0,1.99820031d0,1.01764668d0,
     .	0.32463098d0,2.41041209d0,0.0000005274d0/

	x = R1g-C1(5)
	y = R2g-C2(5)
	ex = exp(-C1(2)*x)
	ey = exp(-C2(2)*y)
	VOH = -C1(1)*(1.0d0 + C1(2)*x + C1(3)*x*x + C1(4)*x*x*x)*ex
	VHCl = -C2(1)*(1.0d0 + C2(2)*y + C2(3)*y*y + C2(4)*y*y*y)*ey
c	VClO = 0.0697d0*exp(-1.1205*(R3g-2.966d0))
	VClO = 0.10302243d0*exp(-0.60616291d0*(R3g-2.966d0))

	return
	end
c ----------------------------------------------------------------------------
	double precision function q(x,xp)
c.. Evaluate the distance-like kernel of order nk with weight function
c.. x**(-mk), given (x,x')

c.. Currently, mk is limited to 6 and nk to 2.
c ----------------------------------------------------------------------------
	implicit double precision (a-h, o-z)
	common/kernel/nk, mk, ns
	dimension F(0:9)

c.. Store a few factorials as floating point numbers.
c.. Calculating these *really* slows things down!
	data (F(i),i=0,9)/1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     .720.0d0, 5040.0d0, 40320.0d0, 362880.0d0/

	if (x .le. xp) then
	    z = x/xp
	    A = 1.0d0/(xp**(mk+1))
        else
	    z = xp/x
	    A = 1.0d0/(x**(mk+1))
        endif

c.. Beta function:
	B = F(nk-1)*F(mk)/F(nk+mk)

c.. The special case of nk = 2 (only case implemented).
	C = float(mk+1)/float(mk+3)
	q = (nk**2)*A*B*(1.0d0 - C*z)

c.. Modification of the distance-like kernel:
c	q = q*(x*xp)**2

	return
     	end
c ----------------------------------------------------------------------------
	double precision function qK(x,xp)
c.. Evaluate the angle-like kernel given (x,x')
c ----------------------------------------------------------------------------
	implicit double precision (a-h, o-z)
	data zero, one/0.0d0, 1.0d0/

	qK = zero
	np = 2
	mp = np-1
	npp = np+1

c.. If this summation is not done, we get the "tilde" kernel.
        do i = 0,mp
	    qK = qK + (x*xp)**i
        enddo

	if (x .le. xp) then
            qK = qK + np*(x**np)*(xp**mp)*(one - mp*x/(npp*xp))
	else
            qK = qK + np*(x**mp)*(xp**np)*(one - mp*xp/(npp*x))
	endif

	return
	end
c-----------------------------------------------------------------------------
	double precision function CFG(a,b,c)

c Given three sides of a triangle A,B,C, find the cosine of the angle A-B-C.
c-----------------------------------------------------------------------------
	implicit double precision (a-h,o-z)

	a2=a**2
	b2=b**2
	c2=c**2
	if (a2.eq.0.0d0) then
	if (b2-c2.lt.0.0d0) then
	x=-1.0d0
	else
	x=1.0d0
	endif
	else if (b2.eq.0.0d0) then
	if (a2-c2.lt.0.0d0) then
	x=-1.0d0
	else
	x=1.0d0
	endif
	else
	x=(a2+b2-c2)/(2.d0*a*b)
	if (x.gt.1.0d0) x=1.0d0
	if (x.lt.-1.0d0) x=-1.0d0
	endif
	CFG=x

	return
	end
c-----------------------------------------------------------------------------
	double precision function cosd(x)

	implicit double precision (a-h,o-z)
	data pi/3.141592653589d0/

	cosd = dcos(x*pi/180.0d0)

	return
	end
