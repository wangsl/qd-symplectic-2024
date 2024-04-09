
/* $Id$ */

/* Reference: 
   Daoqi Yang, C++ and Object-Oriented Numeric Computing for Scientists and Engineers
   Springer-Verlag, New York, Inc. 2001
   P223
*/

#if 0

#ifndef RVECALGO_H
#define RVECALGO_H

class RVec;

struct Scalar_RVec
{
  const double &scalar;
  const RVec &rvec;
  Scalar_RVec(const double &alpha, const RVec &x) : scalar(alpha), rvec(x)
  { };
};

inline Scalar_RVec operator *(const double &alpha, const RVec &x)
{ return Scalar_RVec(alpha, x); }

struct Scalar_RVec_Sum
{
  const Scalar_RVec &ax;
  const Scalar_RVec &by;
  Scalar_RVec_Sum(const Scalar_RVec &ax_, const Scalar_RVec &by_) :
    ax(ax_), by(by_) { }
  Scalar_RVec_Sum(const Scalar_RVec &ax_, const RVec &y_) :
    ax(ax_), by(Scalar_RVec(1.0, y_)) { }
};

inline Scalar_RVec_Sum operator +(const Scalar_RVec &ax, const Scalar_RVec &by)
{ return Scalar_RVec_Sum(ax, by); }

inline Scalar_RVec_Sum operator +(const Scalar_RVec &ax, const RVec &y)
{ return Scalar_RVec_Sum(ax, y); }

inline Scalar_RVec_Sum operator +(const RVec &y, const Scalar_RVec &ax)
{ return Scalar_RVec_Sum(ax, y); }

#endif /* RVECALGO_H */

#endif

