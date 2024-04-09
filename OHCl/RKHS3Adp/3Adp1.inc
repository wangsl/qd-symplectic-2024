
      parameter (n31=255,n32=521)
      common/kernel/nk, mk, ns
      common/geoms/R31(n31,3), R32(n32,3) 
      common/cffts/C31(n31), C32(n32), e0, sf
      data R1a,R2a,thc,dtc/100.0d0,2.41041209d0,75.0d0,16.0d0/
      data zero,one,two,four,eps/0.0d0,1.0d0,2.0d0,4.0d0,1.0d-12/
      data pi/3.14159265358979d0/
      Data e0/-0.1703532D0/

      Character*512 DataDir
      Integer DataDirLength
      Common /RKHSDataDir/DataDirLength, DataDir

      
