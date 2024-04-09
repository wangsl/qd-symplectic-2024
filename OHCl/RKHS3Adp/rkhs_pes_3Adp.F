
C	Modified by Shenglong Wang to run in parallel with OpenMP
C	5/23/2018
	
	subroutine prepot

c   THE 3A" RP SURFACE
c ----------------------------------------------------------------------------
c..					     Version Date: December 6, 2002
C ---DERIVATIVES NOT IMPLEMENTED

c.. The O(3P) + HCl --> OH + Cl MRCI+Q/CBS(aug-cc-pVnZ; n = 2-4) 3A"
c.. potential surface of B. Ramachandran and K. A. Peterson.

c.. The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
c.. to interpolate the surface consisting of 706 geometries spanning
c.. O-H-Cl angles of 51.8-180 deg.  This surface does not contain the
c.. H + ClO arrangement and no claims of accuracy are made for O-H-Cl
c.. angles smaller than 50 deg. 
c ----------------------------------------------------------------------------
c  USAGE:

c  All arguments are passed through the common block /potcm/. The calling
c  routine should have this common block and should pass the array r(3) to
c  the potential subroutine. 
c  On Input:
c	r(1) = R(O-H)
c	r(2) = R(H-Cl)
c	r(3) = R(Cl-O).
c
c On Output:
c	energy = Energy in Hartree w.r to the asymptotic O + HCl minimum.
C ---DERIVATIVES NOT YET IMPLEMENTED
c	dedr   = Derivatives of the potential w.r. to r(1), r(2), and r(3).

c  NOTE:
c    Before any actual potential energy calculations are made, a single
c    call to prepot must be made:
c	call prepot

c    This reads in the three body expansion coefficients, computes the 
c    asymptotic O + HCl minimum energy and saves it. Later, the potential 
c    energy is computed by calling pot:
c	call pot

c ----------------------------------------------------------------------------
	implicit double precision (a-h, o-z)
#if 0	
	parameter (n31=255,n32=521)
	character*40 file1, file2
	common /potcm/r(3), energy, dedr(3)
	common/kernel/nk, mk, ns
	common/geoms/R31(n31,3), R32(n32,3) 
	common/cffts/C31(n31), C32(n32), e0, sf
	common/v2sum/R1g, R2g, R3g, VOH, VHCl, VClO
	data R1a,R2a,thc,dtc/100.0d0,2.41041209d0,75.0d0,16.0d0/
	data zero,one,two,four,eps/0.0d0,1.0d0,2.0d0,4.0d0,1.0d-12/
	data pi/3.141592653589d0/
	external cosd
#endif
	character*40 file1, file2
	Include '3Adp1.inc'

c.. Threebody expansion coefficient files:
	file1 = 'threebody_3Adp_255cm.txt'
	file2 = 'threebody_3Adp_521cm.txt'

	open (2,file=DataDir(1:DataDirLength)//'/'//file1,
	1    form='formatted',status='old')
	read (2,*) sf, nk, mk, ns
	do i = 1,n31
	    read (2,'(3f10.2,d25.15)') (R31(i,j),j=1,3), C31(i)	
	enddo
	close (2)

	open (2,file=DataDir(1:DataDirLength)//'/'//file2,
	1    form='formatted',status='old')
	read (2,*) sf, nk, mk, ns
	do i = 1,n32
	    read (2,'(3f10.2,d25.15)') (R32(i,j),j=1,3), C32(i)	
	enddo
	close (2)

c.. Calculate the asymptotic O + HCl minimum energy e0:
c.. Sum of the two body terms:
#if 0
	R1g = R1a
	R2g = R2a
	R3g = R1a+R2a

	call twobody

	e0 = e0 + VOH + VHCl + VClO
#endif

	write(6,*) 'rkhs_pes_3Adp_n.f'
	write(6,*) 'PREPOT has been called for the OHCl RP 3A" surface.'
	write(6,*) 'Version date December 6, 2002'
	write(6,*) 'Coefficient files:'
	write(6,'(1x,a)') file1
	write(6,'(1x,a)') file2
	write (6,101) nk, mk, ns
	write (6,*) 'Asymptotic Potential: ',e0
	Call Flush(6)
 101	format (1x,
	1    'Order of the kernel n = ',i3/1x,
	2    'Reciprocal power m = ',i3/1x,
	3    'Power of the r coordinate (r**s)  s = ',i3)
 	End

	SubRoutine Pot(rOH, rHCl, rOCl, energy)
	implicit double precision (a-h, o-z)
	Include '3Adp1.inc'
	Include '3Adp2.inc'
	Real*8 r(3)

	r(1) = rOH
	r(2) = rHCl
	r(3) = rOCl

c.. Twobody terms:
	R1g = r(1)
	R2g = r(2)
	R3g = r(3)
	call twobody	

c.. Threebody terms:
	xg = r(1)**ns
	yg = r(2)**ns
	cst = CFG(r(1),r(2),r(3))
	zg = (one - cst)/two + eps
	theta = acos(cst)*180.0d0/pi

	e1 = zero
	e2 = zero

c.. The first threebody term (O-H-Cl angle < 110.0):
	if (theta .le. 110.0d0) then
	    do j = 1,n31
           	xj = R31(j,1)**ns
            	yj = R31(j,2)**ns
	    	zj = (one - cosd(R31(j,3)))/two + eps
           	e1 = e1 + C31(j)*q(xg,xj)*q(yg,yj)*qK(zg,zj)
	    enddo
	endif

c.. The second three body term (O-H-Cl angle > 75):
	if (theta .ge. 75.0d0) then
	    do j = 1,n32
            	xj = R32(j,1)**ns
            	yj = R32(j,2)**ns
	    	zj = (one - cosd(R32(j,3)))/two + eps
            	e2 = e2 + C32(j)*q(xg,xj)*q(yg,yj)*qK(zg,zj)
	    enddo
	endif

c.. The "intermediate" angles:
	e3 = e1 + e2
	if (theta .ge. 75.0d0 .and. theta .le. 110.0d0) then
	    call switch(theta, thc, two, s1)
	    s2 = one - s1
	    e3 = s1*e1 + s2*e2
	endif

c.. The radial switching function:
	dist = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
	if (dist .lt. dtc) then
	    s3 = one
	else
	    call switch (dist, dtc, two, s3)
	endif

c.. Final expression for energy:
	energy = sf*e3*s3 + VOH + VHCl + VClO - e0

	return
#if 0	
101	format (/10x,'Order of the kernel n = ',i3/10x,
     .	'Reciprocal power m = ',i3/10x,
     $'Power of the r coordinate (r**s)  s = ',i3/)
#endif	
	end
c ----------------------------------------------------------------------------
	subroutine twobody

c.. Fit the two body terms using the angle-like kernel and return the value
c.. at R1, R2, R3. On input, R1 = ROH, R2 = RHCl, R3 = RClO in Bohrs.
c -------------------------------------------------------------------------
	implicit double precision (a-h, o-z)
C	common/v2sum/R1g, R2g, R3g, VOH, VHCl, VClO
	Include '3Adp2.inc'
	dimension C1(6), C2(6)
	data (C1(i),i=1,6)/0.17046982d0,2.45278241d0,1.52360833d0,
     .	0.62204843d0,1.83107461d0,0.0000005274d0/
	data (C2(i),i=1,6)/0.17035320d0,1.99820031d0,1.01764668d0,
     .	0.32463098d0,2.41041209d0,0.0000005274d0/

	x = R1g-C1(5)
	y = R2g-C2(5)
	ex = exp(-C1(2)*x)
	ey = exp(-C2(2)*y)
	VOH = -C1(1)*(1.0d0 + C1(2)*x + C1(3)*x*x + C1(4)*x*x*x)*ex
	VHCl = -C2(1)*(1.0d0 + C2(2)*y + C2(3)*y*y + C2(4)*y*y*y)*ey
	VClO = 0.10302243d0*exp(-0.60616291d0*(R3g-2.966d0))

	return
	end
c ----------------------------------------------------------------------------
	double precision function q(x,xp)
c.. Evaluate the distance-like kernel of order nk with weight function
c.. x**(-mk), given (x,x')
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
	C = Dble(mk+1)/Dble(mk+3)
	q = (nk**2)*A*B*(1.0d0 - C*z)

c.. Modification of the distance-like kernel:
	q = q*(x*xp)**2

	return
     	end
c ----------------------------------------------------------------------------
	double precision function qK(x,xp)
c.. Evaluate the angle-like kernel of order np given (x,x')
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
c ----------------------------------------------------------------------------
	double precision function CFG(a,b,c)

c Given three sides of a triangle A,B,C, find the cosine of the angle A-B-C.
c ----------------------------------------------------------------------------
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
c ----------------------------------------------------------------------------
	double precision function cosd(x)

c Calculates the cosine of the angle where the angle is given in degrees
c ----------------------------------------------------------------------------

	implicit double precision (a-h,o-z)
	data pi/3.141592653589d0/

	cosd = dcos(x*pi/180.0d0)

	return
	end
c ----------------------------------------------------------------------------
	subroutine switch (x, x0, fm, s)

c.. The "half-Gaussian" switching function.
c ----------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	data sigma/3.12596d0/

	dx = x - x0
	s = exp(-(dx*dx)/(fm*sigma*sigma))

	return
	end
