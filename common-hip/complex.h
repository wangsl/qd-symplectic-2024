
/* $Id$ */

#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
using namespace std;
#include <cmath>
#include <cstdlib>
#include <cstdio>

#if defined __NVCC__ || defined __HIPCC__
#define __DEVICE_HOST__ __device__ __host__
#else
#define __DEVICE_HOST__
#endif

#ifndef Pi
#define Pi M_PI 
#endif

class Complex
{
 public:
  double re, im;
  
 public:
  
  __DEVICE_HOST__ Complex(double r = 0, double i = 0) : re(r), im(i) { }
  __DEVICE_HOST__ Complex(const Complex &c) : re(c.re), im(c.im) { }
  
  __DEVICE_HOST__ operator double *() { return (double *) this; }
  __DEVICE_HOST__ operator const double *() const { return (const double *) this; }

  __DEVICE_HOST__ double real() const { return re; }
  __DEVICE_HOST__ double imag() const { return im; }
  
  __DEVICE_HOST__ void set_from_polar(double r=0.0, double theta = 0.0)
  { re = r*cos(theta); im = r*sin(theta); }
  
  __DEVICE_HOST__ friend double real(const Complex &c) { return c.re; }
  __DEVICE_HOST__ friend double imag(const Complex &c) { return c.im; }
  __DEVICE_HOST__ friend Complex conj(const Complex &c) { return Complex(c.re, -c.im); }
  
  __DEVICE_HOST__ friend double abs2(const Complex &c) { return c.re*c.re + c.im*c.im; }
  __DEVICE_HOST__ friend double abs(const Complex &c) { return sqrt(abs2(c)); }
  __DEVICE_HOST__ friend double norm(const Complex &c) { return abs(c); }
  
  __DEVICE_HOST__ double norm() const { return abs(*this); }
  __DEVICE_HOST__ double magnitude() const { return abs2(*this); }
  
  __DEVICE_HOST__ friend double arg(const Complex &c)
  {
    double angle = asin(c.im/abs(c));
    if(c.re < 0 && c.im > 0)
      return Pi - angle;
    if(c.re < 0 && c.im < 0)
      return -(Pi+angle);
    return angle;
  }
  
  __DEVICE_HOST__ Complex & operator =(const Complex &c)
  { 
    if(this != &c) {
      re = c.re;
      im = c.im;
    }
    return *this;
  }

  __DEVICE_HOST__ Complex & operator =(const double r)
  { 
    re = r;
    im = 0.0;
    return *this;
  }
  
  __DEVICE_HOST__ friend int operator ==(const Complex &c1, const Complex &c2)
  { return (c1.re == c2.re && c1.im == c2.im); }
  
  __DEVICE_HOST__ friend int operator !=(const Complex &c1, const Complex &c2)
  { return (c1.re != c2.re || c1.im != c2.im); }
  
  __DEVICE_HOST__ Complex operator -() const { return Complex(-re, -im); }
  __DEVICE_HOST__ Complex & operator +=(const Complex &c) { re += c.re; im += c.im; return *this; }
  __DEVICE_HOST__ Complex & operator -=(const Complex &c) { re -= c.re; im -= c.im; return *this; }
  __DEVICE_HOST__ Complex & operator +=(double r) { re += r; return *this; }
  __DEVICE_HOST__ Complex & operator -=(double r) { re -= r; return *this; }
  __DEVICE_HOST__ Complex & operator *=(double r) { re *= r; im *= r; return *this; }
  __DEVICE_HOST__ Complex & operator /=(double r) { return *this *= 1.0/r; }
  __DEVICE_HOST__ Complex & operator *=(const Complex &c) { return *this = *this * c; }
  __DEVICE_HOST__ Complex & operator /=(const Complex &c) { return *this = *this / c; }
  
  __DEVICE_HOST__ friend Complex operator +(double r, const Complex &c) { return Complex(r+c.re, c.im); }
  __DEVICE_HOST__ friend Complex operator +(const Complex &c, double r) { return r + c; }
  __DEVICE_HOST__ friend Complex operator +(const Complex &c1, const Complex &c2)
  { return Complex(c1.re+c2.re, c1.im+c2.im); }
  
  __DEVICE_HOST__ friend Complex operator -(double r, const Complex &c) { return Complex(r-c.re, -c.im); }
  __DEVICE_HOST__ friend Complex operator -(const Complex &c, double r) { return Complex(c.re-r, c.im); }
  __DEVICE_HOST__ friend Complex operator -(const Complex &c1, const Complex &c2)
  { return Complex(c1.re-c2.re, c1.im-c2.im); }
  
  __DEVICE_HOST__ friend Complex operator *(const Complex &c, double r) { return Complex(c.re*r, c.im*r); }
  __DEVICE_HOST__ friend Complex operator *(double r, const Complex &c) { return c*r; }
  __DEVICE_HOST__ friend Complex operator *(const Complex &c1, const Complex &c2)
  { return Complex(c1.re*c2.re - c1.im*c2.im, c1.re*c2.im + c1.im*c2.re); }
  
  __DEVICE_HOST__ friend Complex operator /(const Complex &c, double r) { return Complex(c.re/r, c.im/r); }
  __DEVICE_HOST__ friend Complex operator /(double r, const Complex &c) { return r/abs2(c) * conj(c); }
  __DEVICE_HOST__ friend Complex operator /(const Complex &c1, const Complex &c2) { return c1/abs2(c2) * conj(c2); }
  
  __DEVICE_HOST__ friend Complex exp(const Complex &c) { return exp(c.re) * Complex(cos(c.im), sin(c.im)); }
  __DEVICE_HOST__ friend Complex sinh(const Complex &c) { return (exp(c) - exp(-c)) / 2.0; }
  __DEVICE_HOST__ friend Complex cosh(const Complex &c) { return (exp(c) + exp(-c)) / 2.0; }
  __DEVICE_HOST__ friend Complex tanh(const Complex &c) { return sinh(c)/cosh(c); }
  __DEVICE_HOST__ friend Complex sin(const Complex &c) { return sinh(Complex(0,1.0) * c) / Complex(0,1.0); }
  __DEVICE_HOST__ friend Complex cos(const Complex &c) { return cosh(Complex(0,1.0) * c); }
  __DEVICE_HOST__ friend Complex tan(const Complex &c) { return sin(c)/cos(c); }
  
  __DEVICE_HOST__ friend Complex pow(const Complex &c, double n)
  {
    double theta = n*arg(c);
    double r = pow(abs2(c), n/2.0);
    return Complex(r*cos(theta), r*sin(theta));
  }
  
  __DEVICE_HOST__ friend Complex sqrt(const Complex &c) { return pow(c, 0.5); }
  __DEVICE_HOST__ friend Complex pow(double r, const Complex &c) { return exp(c*log(r)); }

  __DEVICE_HOST__ void zero() { re = 0.0; im = 0.0; }
  
  /********************************************************
   **
   ** Maybe there are something wrong with 'log' operator
   ** when dealing with log(-1) and the similar special
   ** points
   **
   ** *****************************************************/
  
  __DEVICE_HOST__ friend Complex log(const Complex &c)
  {
    double theta = arg(c);
    double rho = abs(c);
    return log(rho) + Complex(0,1.0)*theta;
  }
  
  __DEVICE_HOST__ friend Complex pow(const Complex &c1, const Complex &c2)
  { return exp(c2*log(c1)); }
  
  friend ostream & operator <<(ostream &s, const Complex &c)
  { return s << "(" << c.re << ", " << c.im << ")"; }
  
#if 0
  __DEVICE_HOST__ char *to_string() const 
  { 
    char tmp[32]; 
    sprintf(tmp, "(%.8f, %.8f)", re, im); 
    return tmp;
  }
#endif
};

#endif /* COMPLEX_H */  
