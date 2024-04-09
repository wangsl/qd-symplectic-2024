/* $Id$ */

#ifndef MATBASE_H
#define MATBASE_H

#include <iostream>
using namespace std;
#include <cassert>

template<class T> class Mat
{
  struct MatRep
  {
    int ref, n, m, delete_p;
    T *p;
    T **pp;
    
    MatRep(int nn, int mm) : ref(1), n(nn), m(mm), delete_p(1),
			     p(new T[n*m]), pp(new T*[m])
    { 
      for (int j = 0; j < m; j++)
	pp[j] = p + j*n;
    }
    
    MatRep(int nn, int mm, T *q) : ref(1), n(nn), m(mm), delete_p(0),
				   p(q), pp(new T*[m])
    { 
      for(int j = 0; j < m; j++)
	pp[j] = p + j*n;
    }
    
    ~MatRep() 
    { 
      if(delete_p && p) { delete [] p; p = 0; }
      if(pp) { delete [] pp; pp = 0; }
    } 
    
    T & operator()(int i, int j)
    {
      assert(0 <= i && i < n && 0 <= j && j < m);
      return pp[j][i];
    }
  };
  
  void destroy()
  {
    if (--rep->ref == 0) {
      if(rep) { delete rep; rep = 0; }
    }
  }
  
public:
  MatRep *rep;
  int rows() const
  { return rep->n; }
  int columns() const
  { return rep->m; }
  operator T*()
  { return rep->p; }
  operator const T*() const
  { return rep->p; }
  operator const T*()
  { return rep->p; }
  T & operator()(int i, int j)
  { return (*rep)(i,j); }
  const T & operator()(int i, int j) const
  { return (*rep)(i,j); }
  Mat(int n = 0, int m = 0) : rep(new MatRep(n,m))
  { }
  Mat(const Mat<T> &mat) : rep(mat.rep)
  { rep->ref++; }
  
  Mat(int n, int m, const T &t) : rep(new MatRep(n,m))
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) = t;
  }
  
  Mat(int n, int m, T *q) : rep(new MatRep(n, m, q)) { }
  
  ~Mat()
  { destroy(); }
  
  Mat<T> & operator =(const Mat<T> &mat)
  { 
    mat.rep->ref++;
    destroy();
    rep = mat.rep;
    return *this;
  }
  
  Mat<T> copy() const
  { 
    Mat<T> mat(rows(), columns());
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	mat(i,j) = (*this)(i,j);
    return mat;
  }
  
  void resize(int n, int m)
  {
    if (n == rows() && m == columns())
      return;
    Mat<T> mat(n, m);
    for (int i = 0; i < rows() && i < mat.rows(); i++)
      for (int j = 0; j < columns() && j < mat.columns(); j++)
	mat(i,j) = (*this)(i,j);
    *this = mat;
  }
  
  void resize(int n, int m, const T &t)
  {
    if (n == rows() && m == columns())
      return;
    int i, j;
    Mat<T> mat(n, m);
    for (i = 0; i < mat.rows(); i++)
      for (j = 0; j < mat.columns(); j++)
	if (i < rows() && j < columns())
	  mat(i,j) = (*this)(i,j);
	else
	  mat(i,j) = t;
    *this = mat;
  }
  
  int is_square() const { return rows() == columns(); }
  
  int is_conformable_with(const Mat<T> &mat) const
  { return rows() == mat.rows() && columns() == mat.columns(); }
  
  void symmetrize()
  {
    assert(is_square());
    for (int i = 0; i < rows(); i++)
      for (int j = i+1; j < rows(); j++)
	(*this)(j,i) = (*this)(i,j);
  }
  
  Mat<T> transpose() const
  {
    Mat<T> mat(columns(), rows());
    for (int i = 0; i < rows(); i++) {
      for (int j = 0; j < columns(); j++)
	mat(j,i) = (*this)(i,j);
    }
    return mat;
  }
  
  friend Mat<T> HorizontalPaste(const Mat<T> m1, const Mat<T> m2)
  {
    assert(m1.rows() == m2.rows());
    Mat<T> mat(m1.rows(), m1.columns() + m2.columns());
    for (int i = 0; i < m1.rows(); i++) {
      int j;
      for (j = 0; j < m1.columns(); j++)
	mat(i,j) = m1(i,j);
      for (j = 0; j < m2.columns(); j++)
	mat(i,j+m1.columns()) = m2(i,j);
    }
    return mat;
  }
  
  friend Mat<T> VerticalPaste(const Mat<T> m1, const Mat<T> m2)
  {
    assert(m1.columns() == m2.columns());
    Mat<T> mat(m1.rows() + m2.rows(), m1.columns());
    int i;
    for (i = 0; i < m1.rows(); i++) {
      for (int j = 0; j < m1.columns(); j++)
	mat(i,j) = m1(i,j);
    }
    for (i = 0; i < m2.rows(); i++) {
      for (int j = 0; j < m1.columns(); j++)
	mat(i+m1.rows(), j) = m2(i,j);
    }
    return mat;
  }
};

template<class T> inline ostream & operator<<(ostream &s, const Mat<T> &mat)
{
  s << mat.rows() << " " << mat.columns() << "\n";
  for(int i = 0; i < mat.rows(); i++) {
    for(int j = 0; j < mat.columns(); j++)
      s << mat(i,j) << " ";
    s << "\n";
  }
  return s;
}

template<class T> inline istream & operator>>(istream &s, Mat<T> &mat)
{
  int n, m;
  s >> n >> m;
  if (n != mat.rows() || m != mat.columns())
    mat = Mat<T>(n,m);
  for (int i = 0; i < mat.rows(); i++)
    for (int j = 0; j < mat.columns(); j++)
      s >> mat(i,j);
  return s;
}

#endif /* MATBASE_H */
