
C     Data from Table II in JCP 95, 6421 (1991)
      
      Real*8 a, b, r0, beta0, S, D
      Common /KSGPESData/a(3,3), b(3,3), r0(3), beta0(3), S(3), D(3)
      
      Real*8 r, theta, re, beta
      Common /KSGPESVar/r(3), theta(3), re(3), beta(3)

C     FORT(ohclksgpes) can be run in OpenMP mode only after 
C     declare KSGPESVar as threadprivate

C$OMP THREADPRIVATE(/KSGPESVar/)

