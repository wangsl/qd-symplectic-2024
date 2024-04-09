
/* $Id$ */

#ifndef MATLAB_ARRAY_H
#define MATLAB_ARRAY_H

#include <iostream>
#include <mex.h>
#include "matlabUtils.h"

template<class T> class MatlabArray
{
public:
  MatlabArray(const mxArray *mx_) :
    _mx(mx_), _data(0)
  {    
    _data = (T *) mxGetData(_mx);
    insist(_data);
  }
  
  const mxArray *mx() const { return _mx; }
  const T *data() const { return _data; }
  int n_dims() const  { return _mx ? mxGetNumberOfDimensions(_mx) : 0; }
  const size_t *dims() const { return _mx ? mxGetDimensions(_mx) : 0; }
  
  ~MatlabArray() { _mx = 0; _data = 0; }
  
private:
  const mxArray *_mx;
  T *_data;
};

#ifdef MATLAB2017b
class MatlabComplexArray : public MatlabArray<double>
{
public:

  MatlabComplexArray(const mxArray *mx_) :
    MatlabArray<double>(mx_),
    _imag_data(0)
  {
    insist(mxIsComplex(mx()));
    insist(mxGetPr(mx()) == data());
    
    _imag_data = mxGetPi(mx());
    insist(_imag_data);
  }

  ~MatlabComplexArray() { _imag_data = 0; }

  const double *real() const { return data(); }
  const double *imag() const { return _imag_data; }

private:
  double *_imag_data;
};

#else

// for Matlab after 2018a and after

class MatlabComplexArray
{
public:

  MatlabComplexArray(const mxArray *mx_, const mxArray *my_) :
    x(MatlabArray<double>(mx_)),
    y(MatlabArray<double>(my_))
  {
    insist(mxIsDouble(mx_));
    insist(mxIsDouble(my_));
    check_dims();
  }

  ~MatlabComplexArray() { }

  int n_dims() const  { return x.n_dims(); }
  const size_t *dims() const { return x.dims(); }

  const double *real() const { return x.data(); }
  const double *imag() const { return y.data(); }

private:

  MatlabArray<double> x;
  MatlabArray<double> y;

  void check_dims() const
  {
    insist(x.n_dims() == y.n_dims());

    const size_t *x_dims = x.dims();
    const size_t *y_dims = y.dims();
    for(int i = 0; i < x.n_dims(); i++)
      insist(x_dims[i] == y_dims[i]);
  }
};
#endif

#endif /* MATLABARRAY_H */
