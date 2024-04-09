/* $Id$ */

#ifndef VECBASE_H
#define VECBASE_H

#include "lst.h"
#include "indent.h"

template<class T> class Vec
{
private:
  T *p;
  int n;
  int *ref;
  void destroy()
    { 
      if (ref && --(*ref) == 0) {
	delete[] p;
	delete ref;
	p = 0;
	n = 0;
	ref = 0;
      }
    }
public:
  const int &size() const
    { return n; }
  operator T*()
    { return p; }
  operator const T*() const
    { return p; }
  operator const T*()
    { return p; }
  T & operator[](int i)
    { 
      assert(0 <= i && i < n);
      return p[i];
    }
  const T & operator[](int i) const
    {
      assert(0 <= i && i < n);
      return p[i];
    }
  Vec(int m = 0) : p(new T[m]), n(m), ref(new int(1))
    { }
  Vec(int m, T *q) : p(q), n(m), ref(0)
    { }
  Vec(const Vec<T> &vec) : p(vec.p), n(vec.n), ref(vec.ref)
    { if (ref) (*ref)++; }
  Vec(int m, const T &t) : p(new T[m]), n(m), ref(new int(1))
    {
      for (int i = 0; i < size(); i++)
	(*this)[i] = t;
    }
  Vec(const Vec<T> &v1, const Vec<T> &v2) 
    : p(new T[v1.n+v2.n]), n(v1.n+v2.n), ref(new int(1))
    {
      int i, j;
      for (i = 0; i < v1.size(); i++)
	(*this)[i] = v1[i];
      for (j = 0; j < v2.size(); i++, j++)
	(*this)[i] = v2[j];
    }
  ~Vec()
    { destroy(); }
  Vec<T> & operator =(const Vec<T> &vec)
    { 
      if (vec.ref)
	(*vec.ref)++;
      destroy();
      p = vec.p;
      n = vec.size();
      ref = vec.ref;
      return *this;
    }
  Vec<T> copy() const
    { 
      Vec<T> vec(size());
      for (int i = 0; i < size(); i++)
	vec[i] = (*this)[i];
      return vec;
    }
  Vec<T> & copy(const Vec<T> &vec)
    {
      assert(size() == vec.size());
      for (int i = 0; i < size(); i++)
	(*this)[i] = vec[i];
      return *this;
    }
  void resize(int m)
    {
      if (m == size())
	return;
      Vec<T> vec(m);
      for (int i = 0; i < size() && i < vec.size(); i++)
	vec[i] = (*this)[i];
      *this = vec;
    }
  void resize(int m, const T &t)
    {
      if (m == size())
	return;
      int i;
      Vec<T> vec(m);
      for (i = 0; i < size() && i < vec.size(); i++)
	vec[i] = (*this)[i];
      for ( ; i < vec.size(); i++)
	vec[i] = t;
      *this = vec;
    }
  void add(const T &t)
    {
      resize(size()+1);
      (*this)[size()-1] = t;
    }
  int is_conformable_with(const Vec<T> &v) const
    { return size() == v.size(); }
  void set(const T &t)
    {
      for (int i = 0; i < size(); i++)
	(*this)[i] = t;
    }
  void copy_from_list(const Lst<T> &lst)
    { 
      resize(lst.size());
      LstItr<T> l;
      int i = 0;
      for (l.init(lst); l.ok(); l.next())
	(*this)[i++] = l();
    }

  void zero()
  { memset(p, 0, size()*sizeof(T)); }

  void zeros() { zero(); }

  void show_in_one_line_without_endl() const
  {
    cout << " {";
    for(int i = 0; i < size(); i++)
      cout << " " << (*this)[i];
    cout << " }";
  }

  void show_in_one_line() const
  {
    show_in_one_line_without_endl();
    cout << endl;
  }
};

template<class T> inline ostream & operator<<(ostream &s, const Vec<T> &vec)
{
  s << "{\n";
  IndentPush();
  for (int i = 0; i < vec.size(); i++)
    s << Indent() << vec[i] << "\n";
  IndentPop();
  s << Indent() << "} ";
  return s;
}

template<class T> inline istream & operator>>(istream &s, Vec<T> &vec)
{
  char bracket;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{') {
    Lst<T> l;
    s >> l;
    vec.copy_from_list(l);
  } else {
    int n;
    s >> n;
    if (!s)
      die("Vec: expecting '{' or a number");
    vec.resize(n);
    for (int i = 0; i < n; i++)
      s >> vec[i];
  }
  return s;
}

#endif /* VECBASE_H */
