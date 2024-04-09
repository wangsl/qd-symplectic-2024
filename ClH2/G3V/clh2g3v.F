c
      subroutine prepot
c
c System:          ClH2
c Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
c                  plus E3C (a three-center term)
c                  and modifications to improve the saddle point
c                  bend potential; based on GQQ surface
c
c Interface:       3-1V
c Number of derivatives: 0
c Common name:     G3
c Reference:       Allison, T. C.; Lynch, G. C.; Truhlar, D. G.;
c                  Gordon, M. S.,  J. Phys. Chem., Vol. 100, 13575 (1996).
c
c Note:            This is the modified clh2 potential code.
c                  Vector version, includes no derivatives.
c 
c PREPOT must be called once before any calls to POT are made.  It should
c never be called more than once.
c The potential parameters are included in DATA statements.
c
c The coordinates are defined as follows:
c           R(1) = R(H-Cl)
c           R(2) = R(H-H')
c           R(3) = R(H'-Cl)
c
c The zero of energy for this potential occurs when the Cl atom is
c "infinitely" far from the the H2 diatom, and the H2 diatom is at
c its equilibrium seperation.
c
      implicit double precision (a-h,o-z)

      dimension de(3),re(3),b(3),s(3),pf(3),r(3),y(3),q(3),wj(3),
     &          rr(nt,3),e(nt)

      data de     / 106.447d0, 109.458d0, 106.447d0 /
      data re     / 1.2732d0, 0.74127d0, 1.2732d0 /
      data b      / 1.8674d0, 1.9413d0, 1.8674d0 /
      data s      / 0.1835d0, 0.167d0, 0.1835d0 /
      data rj     / 0.0758016022d0 /
      data ap     / 0.0008627355d0 /
      data at     / 0.2981969860d0 /
      data bp     / 0.1439106543d0 /
      data q2     / 0.6940323070d0 /
      data q4     / 1.6963092005d0 /
      data ckcau  / 627.5095d0 /
      data cangau / 0.529177249d0 /
      data c1     / 3.1053877618397071175758280546d+0 /
      data c2     / -1.8942608491155350536449828878d+0 /
      save

      do 10 i=1,3
        de(i)=de(i)/ckcau
        re(i)=re(i)/cangau
        b(i)=b(i)*cangau
        pf(i)=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
10    continue

      write(6,*) 'PREPOT has been called for the ClH2 surface G3'
      write(6,*) 'Potential last modified on 12 March 1994'

      return

      entry pot(rr,e,nt,nsurf)

      do 30 j=1,nt

        r(1)=rr(j,1)
        r(2)=rr(j,2)
        r(3)=rr(j,3)

        f1=dexp(-ap*(r(1)+r(3))**4)
        f2=dexp(-at*(r(1)-r(3))**2)
        f3=dexp(-bp*(r(1)+r(3)-r(2))**2)
        ctheta=(r(1)*r(1)+r(3)*r(3)-r(2)*r(2))/(2.d0*r(1)*r(3))
        if (ctheta .lt. -1.d0) then
          ctheta=-1.d0
        else if (ctheta .gt. 1.d0) then
          ctheta=1.d0
        end if
        stheta=dsin(dacos(ctheta))
        f4=1.d0+stheta*stheta*(q2+q4*stheta*stheta)
        v3b=f4*(rj*f1*f2*f3)

        x=(r(1)-r(3))/r(2)
        g=(0.5d0*(1.d0+ctheta))**6
        fmod=1.d0+c1*g*(1.d0-x**2)+c2*g*(1.d0-x**4)
        y(1)=dexp(-b(1)*(r(1)-re(1)))
        q(1)=0.5d0*(de(1)*(y(1)-2.d0)+pf(1)*(fmod*y(1)+2.d0))*y(1)
        wj(1)=0.5d0*(de(1)*(y(1)-2.d0)-pf(1)*(fmod*y(1)+2.d0))*y(1)
        y(2)=dexp(-b(2)*(r(2)-re(2)))
        q(2)=0.5d0*(de(2)*(y(2)-2.d0)+pf(2)*(y(2)+2.d0))*y(2)
        wj(2)=0.5d0*(de(2)*(y(2)-2.d0)-pf(2)*(y(2)+2.d0))*y(2)
        y(3)=dexp(-b(3)*(r(3)-re(3)))
        q(3)=0.5d0*(de(3)*(y(3)-2.d0)+pf(3)*(fmod*y(3)+2.d0))*y(3)
        wj(3)=0.5d0*(de(3)*(y(3)-2.d0)-pf(3)*(fmod*y(3)+2.d0))*y(3)

        rad=dsqrt(wj(1)*(wj(1)-wj(2))+wj(2)*(wj(2)-wj(3))
     &            +wj(3)*(wj(3)-wj(1)))
        e(j)=de(2)+q(1)+q(2)+q(3)-rad+v3b

30    continue

      return

      end
