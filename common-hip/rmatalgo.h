
/* $Id$ */

#if 0

#ifndef RMATALOG_H
#define RMATALOG_H

#include "rvecalgo.h"

class RVec;
class RMat;

struct Scalar_RMat
{
  const double &scalar;
  const RMat &rmat;
  Scalar_RMat(const double &alpha, const RMat &A) : scalar(alpha), rmat(A) { }
};

inline Scalar_RMat operator *(const double &alpha, const RMat &A)
{ return Scalar_RMat(alpha, A); }

struct Scalar_RMat_RVec
{
  const Scalar_RMat &scalar_rmat;
  const RVec &rvec;
  Scalar_RMat_RVec(const Scalar_RMat &aA, const RVec &x) :
    scalar_rmat(aA), rvec(x) { }
};

inline Scalar_RMat_RVec operator *(const Scalar_RMat &aA, const RVec &x)
{ return Scalar_RMat_RVec(aA, x); }

inline Scalar_RMat_RVec operator *(const RMat &A, const RVec &x)
{ return Scalar_RMat_RVec(Scalar_RMat(1.0, A), x); }

struct Scalar_RVec_RMat
{
  const Scalar_RVec &scalar_rvec;
  const RMat &rmat;
  Scalar_RVec_RMat(const Scalar_RVec &ax, const RMat &A) : 
    scalar_rvec(ax), rmat(A) { }
};

inline Scalar_RVec_RMat operator *(const Scalar_RVec &ax, const RMat &A)
{ return Scalar_RVec_RMat(ax, A); }

inline Scalar_RVec_RMat operator *(const RVec &x, const RMat &A)
{ return Scalar_RVec_RMat(Scalar_RVec(1.0, x), A); }

#endif /* RMATALOG_H */

#endif

