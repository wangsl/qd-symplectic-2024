

C     the subroutine gives the analytical expression for
C     O(3P)+HCl(X 1Sigma) --> OH(X 2Pi)+Cl(2P)

C     H. Koizumi, George C. Schatz, Mark S. Gordon
C     J.Chem.Phys. 95, 6421, (1991)

C     Lichang Wang, C.Kalyanaraman, Anne B.McCoy
C     J.Chem.Phys. 110, 11221, (1999)

      Function KSG_S1(a, x, y)
      Implicit None
      Real*8 KSG_S1, a, x, y
      Real*8 S1
      S1 = 0.5D0 * (1.0D0 + Tanh(a*(x-y)))
      KSG_S1 = S1
      End
      
      Function KSG_S2(a, x)
      Implicit None
      Real*8 KSG_S2, a, x
      Real*8 S2
      S2 = Exp(-a*x*x)
      KSG_S2 = S2
      End
      
      Function KSG_re(i)
      Implicit None
      Real*8 KSG_re
      Integer i
      Include 'OHClKSG.inc'
      Real*8 theta_i, theta_i_2, theta_i_4
      Real*8 S1, S2, re_i
      Real*8 KSG_S1, KSG_S2
      External KSG_S1, KSG_S2
      
      theta_i = theta(i)
      theta_i_2 = theta_i*theta_i
      theta_i_4 = theta_i_2*theta_i_2

      S1 = KSG_S1(2.0D0, 9.0D0, r(2))
      S2 = KSG_S2(6.0D0, Abs(theta_i))
      
      re_i = r0(i) + 
     $     (a(i,1) + a(i,2)*theta_i_2 + a(i,3)*theta_i_4 - r0(i))
     $     * S2 * S1
      
      KSG_re = re_i
      End

      Function KSG_beta(i)
      Implicit None
      Integer I
      Real*8 KSG_beta
      Include 'OHClKSG.inc'
      Real*8 theta_i, theta_i_2, theta_i_4
      Real*8 S1, S2, beta_i

      Real*8 KSG_S1, KSG_S2
      External KSG_S1, KSG_S2
      
      theta_i = theta(i)
      theta_i_2 = theta_i*theta_i
      theta_i_4 = theta_i_2*theta_i_2

      S1 = KSG_S1(2.0D0, 9.0D0, r(2))
      S2 = KSG_S2(3.0D0, Abs(theta_i))

      beta_i = beta0(i) +
     $     (b(i,1) + b(i,2)*theta_i_2 +
     $     b(i,3)*theta_i_4 - beta0(i))
     $     * S2 * S1
      
      KSG_beta = beta_i 
      End

      Function KSG_Q(i)
      Implicit None
      Real*8 KSG_Q
      Integer i
      Include 'OHClKSG.inc'

      Real*8 e_b_r_ri, beta_i, Q_i

      beta_i = beta(i)

      e_b_r_ri = Exp(-beta_i * (r(i)-re(i)))
      
      Q_i = 0.25D0*D(i)*( 
     $     (3.0D0+S(i))*e_b_r_ri*e_b_r_ri 
     $     - (2.0D0+6.0D0*S(i))*e_b_r_ri)/(1.0D0+S(i))

      KSG_Q = Q_i
      End

      Function KSG_J(i)
      Implicit None
      Real*8 KSG_J
      Integer i
      Include 'OHClKSG.inc'
      Real*8 e_b_r_ri, beta_i, J_i
      
      beta_i = beta(i)

      e_b_r_ri = Exp(-beta_i*(r(i) - re(i)))
      
      J_i = 0.25D0*D(i)*( (1.0D0+3.0D0*S(i))*e_b_r_ri*e_b_r_ri -
     $     (6.0D0+2.0D0*S(i))*e_b_r_ri)/(1.0D0+S(i))

      KSG_J = J_i
      End

      
      SubRoutine OHClKSGPES(rOH, rOCl, rHCl, V)
      Implicit None
      Real*8 rOH, rOCl, rHCl, V
      Include 'OHClKSG.inc'
      Real*8 Q1, Q2, Q3, J1, J2, J3
      Real*8 KSG_re, KSG_beta, KSG_Q, KSG_J
      External KSG_re, KSG_beta, KSG_Q, KSG_J

      Call InitKSHParameters

      r(1) = rOH
      r(2) = rOCl
      r(3) = rHCl

      theta(1) = ACos((r(1)*r(1) + r(2)*r(2) - r(3)*r(3))/(2*r(1)*r(2)))
      theta(3) = ACos((r(2)*r(2) + r(3)*r(3) - r(1)*r(1))/(2*r(2)*r(3)))
      theta(2) = theta(1) + theta(3)

      re(1) = KSG_re(1)
      re(2) = KSG_re(2)
      re(3) = KSG_re(3)

      beta(1) = KSG_beta(1)
      beta(2) = KSG_beta(2)
      beta(3) = KSG_beta(3)

      Q1 = KSG_Q(1)
      Q2 = KSG_Q(2)
      Q3 = KSG_Q(3)

      J1 = KSG_J(1)
      J2 = KSG_J(2)
      J3 = KSG_J(3)

      V = Q1 + Q2 + Q3 -
     $     Sqrt(J1*J1 + J2*J2 + J3*J3 -
     $     J1*J2 - J2*J3 - J3*J1)

      End
      
      SubRoutine InitKSHParameters
      Implicit None
      Include 'OHClKSG.inc'
      
      Logical HasParameters
      Data HasParameters /.False./
      Save HasParameters

      If(HasParameters) Return

C     Data from Table II in JCP 95, 6421 (1991)

      a(1,1) = 1.705D0
      a(2,1) = 3.273D0
      a(3,1) = 2.306D0
      a(1,2) = 0.7823D0 
      a(2,2) = -0.9621D0 
      a(3,2) = 0.4932D0
      a(1,3) = 0.4463D0
      a(2,3) = -3.116D0
      a(3,3) = 2.131D0
      
      b(1,1) = 1.313D0
      b(2,1) = 1.268D0
      b(3,1) = 1.239D0
      b(1,2) = -1.206D0
      b(2,2) = -2.934D0
      b(3,2) = -2.840D0
      b(1,3) = 1.763D0
      b(2,3) = 4.700D0
      b(3,3) = 6.012D0

      r0(1) = 1.833D0
      r0(2) = 2.967D0
      r0(3) = 2.409D0
      
      beta0(1) = 1.214D0
      beta0(2) = 1.212D0
      beta0(3) = 0.9885D0
      
      S(1) = 0.0355D0
      S(2) = 0.307D0
      S(3) = 0.153D0

      D(1) = 0.1699D0
      D(2) = 0.1030D0
      D(3) = 0.1697D0

      End
      
      SubRoutine OHClKSGPES2(R, V)
      Implicit None
      Real*8 R(3), V
      Call OHClKSGPES(R(1), R(2), R(3), V)
      End


