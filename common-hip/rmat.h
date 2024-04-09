/* $Id$ */

#ifndef RMAT_H
#define RMAT_H

#include <string.h>
#include <cassert>
#include "matbase.h"
#include "vecbase.h"
#include "fort.h"
#include "rvecalgo.h"
#include "rmatalgo.h"
#include "fint.h"

extern "C" 
{
  void FORT(rvadd)(double*, const double*, const double*, const FInt &);
  void FORT(rvsubt)(double*, const double*, const double*, const FInt &);
  void FORT(rsvmult)(double*, const double*, const double*, const FInt &);
  void FORT(rvvmult)(double*, const double*, const double*, const FInt &);
  void FORT(rmadd)(double*, const double*, const double*, const FInt &, const FInt &);
  void FORT(rmsubt)(double*, const double*, const double*, const FInt &, const FInt &);
  void FORT(rsmmult)(double*, const double*, const double*, const FInt &, const FInt &);
  void FORT(rvmmult)(double*, const double*, const double*, const FInt &, const FInt &);
  void FORT(rmvmult)(double*, const double*, const double*, const FInt &, const FInt &);
  void FORT(rmmmult)(double*, const double*, const double*, const FInt &, const FInt &,
		     const FInt &);
}

class RMat;

class RVec : public Vec<double>
{
public:
  RVec(int sz = 0) : Vec<double>(sz) { } 
  RVec(const RVec &v) : Vec<double>(v) { }
  RVec(const Vec<double> &v) : Vec<double>(v) { }
  RVec(int sz, const double &t) : Vec<double>(sz,t) { }

  RVec(int sz, double *q) : Vec<double>(sz, q) { }
  RVec(double *q, int sz) : Vec<double>(sz, q) { }

  //By Shenglong Wang 2004-01-25
#if 0
  RVec(const Scalar_RVec &ax);
  RVec & operator =(const Scalar_RVec &ax);
  RVec(const Scalar_RVec_Sum &axby);
  RVec & operator =(const Scalar_RVec_Sum &axby);
  RVec(const Scalar_RVec_RMat &axA);
  RVec & operator =(const Scalar_RVec_RMat &axA);
  RVec(const Scalar_RMat_RVec &aAx);
  RVec & operator =(const Scalar_RMat_RVec &aAx);
#endif
  
  //end 2004-01-25
  
  RVec & operator=(const RVec &v)
    { 
      Vec<double>::operator=(v);
      return *this;
    }
  RVec & apply(double (*f)(double))
    { 
      for (int i = 0; i < size(); i++)
	(*this)[i] = f((*this)[i]);
      return *this;
    }
  RVec map(double (*f)(double)) const
    {
      RVec vec(size());
      for (int i = 0; i < size(); i++)
	vec[i] = f((*this)[i]);
      return vec;
    }
  RVec & add(const RVec &a, const RVec &b)
    {
      assert(is_conformable_with(a));
      assert(a.is_conformable_with(b));
      FORT(rvadd)(*this, a, b, size());
      return *this;
    }
  RVec & subtract(const RVec &a, const RVec &b)
    {
      assert(is_conformable_with(a));
      assert(a.is_conformable_with(b));
      FORT(rvsubt)(*this, a, b, size());
      return *this;
    }
  RVec & multiply(double a, const RVec &b)
    {
      assert(is_conformable_with(b));
      FORT(rsvmult)(*this, &a, b, size());
      return *this;
    }
  RVec & multiply(const RVec &, const RMat &);
  RVec & multiply(const RMat &, const RVec &);
  RVec & symmetric_multiply(const RMat &, const RVec &);
  RVec & operator +=(const RVec &t)
  { return add((const RVec &) *this, t); }
  RVec & operator -=(const RVec &t)
  { return subtract((const RVec &) *this, t); }
  RVec & operator *=(double r)
  { return multiply(r, (const RVec &) *this); }
  RVec & operator /=(double r)
  { return multiply(1.0/r, (const RVec &) *this); }
  RVec operator -() const
  { 
    RVec v(size());
    return v.multiply(-1.0, *this);
  }
  RVec & copy(const RVec &vec)
  {
    assert(is_conformable_with(vec));
    memcpy((double *) *this, (const double *) vec, size()*sizeof(double));
    return *this;
  }
  RVec & copy(const double *vec)
  { 
    memcpy((double *) *this, vec, size()*sizeof(double));
    return *this;
  }
  const RVec & copy_to(double *vec) const
  { 
    memcpy(vec, (const double *) *this, size()*sizeof(double));
    return *this;
  }
  RVec & copy(const RVec &vec, int start, int end)
  {
    assert(size() >= end - start + 1);
    assert(0 <= start && start < vec.size());
    assert(0 <= end && end < vec.size());
    memcpy(*this, &vec[start], (end - start + 1) * sizeof(double));
    return *this;
  }
  RVec copy() const
  {
    RVec r(size());
    r.copy(*this);
    return r;
  }
  /*
    void zero()
    { memset(*this, 0, size()*sizeof(double)); }
  */
  
  void garbage() { memset(*this, 0xFF, size()*sizeof(double)); }
  
  double sum() const
  {
    double s = 0.0;
    for (int i = 0; i < size(); i++)
      s += (*this)[i];
    return s;
  }
  
  void show_self();
  
  double sq() const
  {
    double c;
    FORT(rvvmult)(&c, *this, *this, size());
    return c;
  }
  
  double minimum() const
  {
    if(size() == 0) return 0;
    double m = (*this)[0];
    for (int i = 1; i < size(); i++)
      if ((*this)[i] < m)
	m = (*this)[i];
    return m;
  }
  
  double maximum() const
  {
    if (size() == 0) return 0;
    double m = (*this)[0];
    for (int i = 1; i < size(); i++)
      if ((*this)[i] > m)
	m = (*this)[i];
    return m;
  }
  
  void append(const double c) 
  { Vec<double>::add(c); }
  
};

inline RVec operator+(const RVec &a, const RVec &b)
{
  RVec c(a.size());
  return c.add(a,b);
}

inline RVec operator-(const RVec &a, const RVec &b)
{
  RVec c(a.size());
  return c.subtract(a,b);
}

/*
  inline RVec operator*(double a, const RVec &b)
  {
  RVec c(b.size());
  return c.multiply(a,b);
  }
  
  inline RVec operator*(const RVec &a, double b)
  {
  return b*a;
  }
*/

inline double operator *(const RVec &a, const RVec &b)
{
  double c;
  assert(a.is_conformable_with(b));
  FORT(rvvmult)(&c, a, b, a.size());
  return c;
}

class RMat : public Mat<double>
{
public:
  RMat(int n = 0, int m = 0) : Mat<double>(n, m) { }
  RMat(const RMat &mat) : Mat<double>(mat) { }
  RMat(const Mat<double> &mat) : Mat<double>(mat) { }
  RMat(int n, int m, const double &t) : Mat<double>(n, m, t) { }
  RMat(int n, int m, double *q) : Mat<double>(n, m, q) { }

  //By Shenglong Wang 2004-01-25
#if 0
  RMat(const Scalar_RMat &aA);  
  RMat & operator =(const Scalar_RMat &aA);
#endif
  //end 2004-01-25

  RMat & operator=(const RMat &mat)
    { 
      Mat<double>::operator=(mat); 
      return *this;
    }
  RMat transpose() const 
    { return RMat(Mat<double>::transpose()); }
  RMat & copy(const RMat &mat)
    {
      assert(is_conformable_with(mat));
      memcpy(rep->p, mat.rep->p, rows()*columns()*sizeof(double));
      return *this;
    }
  RMat copy() const
    {
      RMat r(rows(), columns());
      r.copy(*this);
      return r;
    }
  int is_symmetric() const
    {
      if (!is_square())
	return 0;
      for (int i = 0; i < rows(); i++)
	for (int j = i+1; j < columns(); j++)
	  if ((*this)(i,j) != (*this)(j,i))
	    return 0;
      return 1;
    }
  double minimum() const
  {
    if (rows() == 0 || columns() == 0)
      return 0;
    double min = (*this)(0,0);
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	if ((*this)(i,j) < min)
	  min = (*this)(i,j);
    return min;
  }
  double maximum() const
  {
    if (rows() == 0 || columns() == 0)
      return 0;
    double max = (*this)(0,0);
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	if ((*this)(i,j) > max)
	  max = (*this)(i,j);
    return max;
  }
  double minimum_in_row(int i) const
  {
    assert(0 <= i && i < rows());
    if (columns() == 0)
      return 0;
    double min = (*this)(i,0);
    for (int j = 1; j < columns(); j++)
      if ((*this)(i,j) < min)
	min = (*this)(i,j);
    return min;
  }
  double maximum_in_row(int i) const
  {
    assert(0 <= i && i < rows());
    if (columns() == 0)
      return 0;
    double max = (*this)(i,0);
    for (int j = 1; j < columns(); j++)
      if ((*this)(i,j) > max)
	max = (*this)(i,j);
    return max;
  }
  double minimum_in_column(int j) const
  {
    assert(0 <= j && j < columns());
    if (rows() == 0)
      return 0;
    double min = (*this)(0,j);
    for (int i = 1; i < rows(); i++)
      if ((*this)(i,j) < min)
	min = (*this)(i,j);
    return min;
  }
  double maximum_in_column(int j) const
  {
    assert(0 <= j && j < columns());
    if (rows() == 0)
      return 0;
    double max = (*this)(0,j);
    for (int i = 1; i < rows(); i++)
      if ((*this)(i,j) > max)
	max = (*this)(i,j);
    return max;
  }
  RMat & apply(double (*f)(double))
    {
      for (int i = 0; i < rows(); i++)
	for (int j = 0; j < columns(); j++)
	  (*this)(i,j) = f((*this)(i,j));
      return *this;
    }
  RMat & add(const RMat &a, const RMat &b)
    {
      assert(is_conformable_with(a));
      assert(a.is_conformable_with(b));
      FORT(rmadd)(*this, a, b, rep->n, rep->m);
      return *this;
    }
  RMat & subtract(const RMat &a, const RMat &b)
    {
      assert(is_conformable_with(a));
      assert(a.is_conformable_with(b));
      FORT(rmsubt)(*this, a, b, rep->n, rep->m);
      return *this;
    }
  RMat & multiply(double a, const RMat &b)
    {
      assert(is_conformable_with(b));
      FORT(rsmmult)(*this, &a, b, rep->n, rep->m);
      return *this;
    }
  RMat & multiply(const RMat &a, const RMat &b)
    {
      assert(a.rep != rep && b.rep != rep);
      assert(rows() == a.rows() && columns() == b.columns());
      assert(a.columns() == b.rows());
      FORT(rmmmult)(*this, a, b, rep->n, rep->m, a.rep->m);
      return *this;
    }
  RMat & l_symmetric_multiply(const RMat &a, const RMat &b);
  RMat & operator+=(const RMat &t)
    { return add((const RMat &) *this, t); }
  RMat & operator-=(const RMat &t)
    { return subtract((const RMat &) *this, t); }
  RMat & operator*=(double r)
    { return multiply(r, (const RMat &) *this); }
  RMat operator-() const
    { 
      RMat m(rows(), columns());
      return m.multiply(-1.0, *this);
    }
  void zero()
    { memset(rep->p, 0, rows()*columns()*sizeof(double)); }
  void garbage()
    { memset(rep->p, 0xFF, rows()*columns()*sizeof(double)); }
  void show_self();
};

inline RMat operator+(const RMat &a, const RMat &b)
{
  RMat c(a.rows(), a.columns());
  return c.add(a,b);
}

inline RMat operator-(const RMat &a, const RMat &b)
{
  RMat c(a.rows(), a.columns());
  return c.subtract(a,b);
}

/*
  inline RMat operator*(double a, const RMat &b)
  {
  RMat c(b.rows(), b.columns());
  return c.multiply(a,b);
  }
  
  inline RMat operator*(const RMat &a, double b)
  {
  return b*a;
  }
*/

inline RMat operator*(const RMat &a, const RMat &b)
{
  RMat c(a.rows(), b.columns());
  return c.multiply(a,b);
}

/*
inline RVec operator*(const RVec &a, const RMat &b)
{
  RVec c(b.columns());
  return c.multiply(a,b);
}

inline RVec operator*(const RMat &a, const RVec &b)
{
  RVec c(a.rows());
  return c.multiply(a,b);
}
*/

inline ostream & operator <<(ostream &s, const RVec &v)
{
  s << v.size();
  if(v.size()) s << "\n";
  for (int i = 0; i < v.size(); i++) {
    s << " " << v[i];
    if(i != v.size()-1) s << "\n";
  }
  return s;
}

inline ostream & operator <<(ostream &s, const RMat &m)
{
  s << m.rows() << " " << m.columns() << "\n";
  for (int i = 0; i < m.rows(); i++) {
    s << " ";
    for (int j = 0; j < m.columns(); j++)
      s << m(i,j) << " ";
    s << "\n";
  }
  return s;
}

#endif /* RMAT_H */
