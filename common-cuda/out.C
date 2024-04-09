/* $Id$ */

#include <iostream>
using namespace std;
#include <fstream>
#include "out.h"
#include "die.h"
#include "fort.h"

static ostream *out = 0;
static int precision = 4;

void OutInit() 
{ 
  SetPrecision(precision);
}

void OutSetToFixed()
{
  Out().setf(ios::fixed);
}

void OutUnsetToFixed()
{
  Out().unsetf(ios::fixed);
}

void OutSetToFile(const char *f) 
{
  if (out)
    delete out;
  out = FileStream(f);
}

void SetPrecision(int p)
{
  Out().precision(precision = p);
}

ostream &Out()
{
  return out ? *out : cout;
}

ostream *FileStream(const char *fname)
{
  ostream *f = new ofstream(fname);
  if (!(f && *f)) {
    Out() << "FileStream: cannot create " << fname << "\n" << flush;
    die("");
  }
  f->precision(precision);
  return f;
}

istream *InFileStream(const char *fname)
{
  istream *f = new ifstream(fname);
  if (!(f && *f)) {
    Out() << "InFileStream: cannot open " << fname << "\n" << flush;
    die("");
  }
  return f;
}



