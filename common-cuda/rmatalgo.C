
/* $Id$ */

#if 0

#include "rmat.h"

RVec::RVec(const Scalar_RVec &ax) : Vec<double>(ax.rvec.size())
{
  multiply(ax.scalar, ax.rvec);
}

RVec & RVec::operator =(const Scalar_RVec &ax)
{
  if(size() != ax.rvec.size())
    resize(ax.rvec.size());
  multiply(ax.scalar, ax.rvec);
  return *this;
}

RVec::RVec(const Scalar_RVec_Sum &axby) : Vec<double>(axby.ax.rvec.size())
{
  assert(axby.ax.rvec.size() == axby.by.rvec.size());
  for(int i = 0; i < size(); i++)
    (*this)[i] = axby.ax.scalar*axby.ax.rvec[i] + axby.by.scalar*axby.by.rvec[i];
}

RVec & RVec::operator =(const Scalar_RVec_Sum &axby)
{
  assert(axby.ax.rvec.size() == axby.by.rvec.size());
  if(size() != axby.ax.rvec.size())
    resize(axby.ax.rvec.size());
  for(int i = 0; i < size(); i++)
    (*this)[i] = axby.ax.scalar*axby.ax.rvec[i] + axby.by.scalar*axby.by.rvec[i];
  return *this;
}

RMat::RMat(const Scalar_RMat &aA) : Mat<double>(aA.rmat.rows(), aA.rmat.columns())
{
  multiply(aA.scalar, aA.rmat);
}

RMat & RMat::operator =(const Scalar_RMat &aA)
{
  if(!is_conformable_with(aA.rmat))
    resize(aA.rmat.rows(), aA.rmat.columns());
  multiply(aA.scalar, aA.rmat);
  return *this;
}

RVec::RVec(const Scalar_RVec_RMat &axA) : Vec<double>(axA.rmat.columns())
{
  assert(axA.scalar_rvec.rvec.size() == axA.rmat.rows());
  const int n = axA.rmat.rows();
  const int m = axA.rmat.columns();
  for(int i = 0; i < m; i++) {
    (*this)[i] = 0.0;
    for(int j = 0; j < n; j++)
      (*this)[i] += axA.scalar_rvec.rvec[j]*axA.rmat(j,i);
    (*this)[i] *= axA.scalar_rvec.scalar;
  }
}

RVec & RVec::operator =(const Scalar_RVec_RMat &axA)
{
  assert(size() == axA.rmat.columns());
  if(size() != axA.rmat.columns())
    resize(axA.rmat.columns());
  const int n = axA.rmat.rows();
  const int m = axA.rmat.columns();
  for(int i = 0; i < m; i++) {
    (*this)[i] = 0.0;
    for(int j = 0; j < n; j++)
      (*this)[i] += axA.scalar_rvec.rvec[j]*axA.rmat(j,i);
    (*this)[i] *= axA.scalar_rvec.scalar;
  }
  return *this;
}

RVec::RVec(const Scalar_RMat_RVec &aAx) : Vec<double>(aAx.scalar_rmat.rmat.rows())
{
  assert(aAx.scalar_rmat.rmat.columns() == aAx.rvec.size());
  const int n = aAx.scalar_rmat.rmat.rows();
  const int m = aAx.scalar_rmat.rmat.columns();
  for(int i = 0; i < n; i++) {
    (*this)[i] = 0.0;
    for(int j = 0; j < m; j++)
      (*this)[i] += aAx.scalar_rmat.rmat(i,j)*aAx.rvec[j];
    (*this)[i] *= aAx.scalar_rmat.scalar;
  }
}

RVec & RVec::operator =(const Scalar_RMat_RVec &aAx)
{
  assert(aAx.scalar_rmat.rmat.columns() == aAx.rvec.size());
  if(size() != aAx.scalar_rmat.rmat.rows())
    resize(aAx.scalar_rmat.rmat.rows());
  const int n = aAx.scalar_rmat.rmat.rows();
  const int m = aAx.scalar_rmat.rmat.columns();
  for(int i = 0; i < n; i++) {
    (*this)[i] = 0.0;
    for(int j = 0; j < m; j++)
      (*this)[i] += aAx.scalar_rmat.rmat(i,j)*aAx.rvec[j];
    (*this)[i] *= aAx.scalar_rmat.scalar;
  }
  return *this;
}

#endif

