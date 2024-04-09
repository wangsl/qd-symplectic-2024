
/* $Id$ */

/**********************************************
  Note by Shenglong Wang, 2-13-2006
  the Integer variables for 64-bits compiler 
  -I64 or -qintsize=8

  IBM compiler: xlc_r xlcc_r xlf_r

  Fortran sizef(integer) = 8
  Fortran: sizeof(Logical) = 8
  Fortran: sizeof(Real*8) = 8
  Fortran: sizeof(Double Precision) = 16

  C++: sizeof(long) = 8
  C++: sizeof(int) = 4
  C++: sizeof((void *) 0) = 8, pointer

**********************************************/

#ifndef FINT_H
#define FINT_H

#ifdef I64
typedef long int FInt;
typedef long int INTEGER;

#else

typedef int FInt;
typedef int INTEGER;

#endif /* I64 */

#endif /* FINT_H */
