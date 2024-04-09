
/* $Id$ */

#include <iostream>
using namespace std;
#include <cassert>
#include <cstring>
#include <cassert>
#include <mex.h>
#include <string.h>
#include "matlabUtils.h"

int file_exist(const string file_name)
{
  return access(file_name.c_str(), F_OK) ? 0 : 1;
}

void mex_to_matlab(const char *script, const int nrhs, mxArray *prhs[])
{
  if(!file_exist(script + string(".m"))) return;

  cout << " Matlab script " << script << endl;

  insist(!mexCallMATLAB(0, NULL, nrhs, prhs, script));
}

void mex_to_matlab(const char *script)
{
  if(!file_exist(script + string(".m"))) return;

  cout << " Matlab script " << script << endl;

  insist(!mexCallMATLAB(0, NULL, 0, NULL, script));
}
