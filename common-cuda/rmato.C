/* $Id$ */

#include "rmat.h"
#include "out.h"
#include "fort.h"

void RVec::show_self() { Out() << *this << "\n" << flush; }
void RMat::show_self() { Out() << *this << "\n" << flush; }

static const char jobL[128] = "L";
static const char jobU[128] = "U";

extern "C" void FORT(dsymm)(const char *, const char *, const FInt *, const FInt *, 
			    const double *, const double *, const FInt *, const double *, 
			    const FInt *, const double *, double *, const FInt *);
extern "C" void FORT(dsymv)(const char *, const FInt &, const double &, const double *, 
			    const FInt &, const double *, const FInt &, const double &, 
			    double *, const FInt &);

RMat & RMat::l_symmetric_multiply(const RMat &a, const RMat &b)
{
  double one(1), zero(0);
  FInt m = b.rows();
  FInt n = b.columns();
  FInt ld = m > 0 ? m : 1;
  assert(a.rows() == m && a.columns() == m);
  assert(rows() == m && columns() == n);
  assert((const double *) a != (const double *) (*this));
  assert((const double *) b != (const double *) (*this));
  FORT(dsymm)(jobL, jobU, &m, &n, &one, a, &ld, b, &ld, &zero, *this, &ld);
  return *this;
}

RVec & RVec::symmetric_multiply(const RMat &a, const RVec &b)
{
  assert(size() == a.rows());
  assert(is_conformable_with(b));
  assert(a.is_square());
  FInt lda = size() > 1 ? size() : 1;
  assert((const double *) b != (const double *) (*this));
  FInt size_l = size();
  FORT(dsymv)(jobU, size_l, 1.0, a, lda, b, 1, 0.0, *this, 1);
  return *this;
}

RVec & RVec::multiply(const RVec &a, const RMat &b)
{
  assert(size() == b.columns());
  assert(a.size() == b.rows());
  assert((const double *) a != (const double *) (*this));
  FInt n = b.rep->n;
  FInt m = b.rep->m;
  FORT(rvmmult)(*this, a, b, n, m);
  return *this;
}

RVec & RVec::multiply(const RMat &a, const RVec &b)
{
  assert(size() == a.rows());
  assert(b.size() == a.columns());
  assert((const double *) b != (const double *) (*this));
  FInt n = a.rep->n;
  FInt m = a.rep->m;
  FORT(rmvmult)(*this, a, b, n, m);
  return *this;
}
