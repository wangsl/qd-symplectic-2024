
/* $Id$ */

#ifndef MATLAB_UTILS_H
#define MATLAB_UTILS_H

#define MCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatlabCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)

/***
 * do not use printf macro in matlab 
 * printf is defined as mexPrintf in Matlab mex.h
 ***/
#ifdef printf
#undef printf
#endif

#define _mxInt32_ 32
#define _mxDouble_ 64

#ifdef __cplusplus

#include <mex.h>
#include <cstring>
#include <unistd.h>

#define insist(x) if (!(x)) MatlabCrashLoc("insist failed: " #x, __FILE__, __LINE__)

void MatlabCrashLoc(const char *message, const char *file_name, const int line);

void mex_to_matlab(const char *script, const int nrhs, mxArray *prhs[]);
void mex_to_matlab(const char *script);

inline int file_exist(const char *file_name)
{ return access(file_name, F_OK) ? 0 : 1; }

inline void *mxGetData(const mxArray *mx, const char *field, const int type = _mxDouble_)
{
  insist(mx);
  insist(type == _mxInt32_ || type == _mxDouble_);
  
  mxArray *mxPtr = mxGetField(mx, 0, field);
  if(!mxPtr) return 0;
  
  if(type == _mxInt32_) {
    insist(mxIsInt32(mxPtr));
  } else if(type == _mxDouble_) {
    insist(mxIsDouble(mxPtr));
  }    
  
  void *ptr = mxGetData(mxPtr);
  return ptr;
}

inline char *mxGetString(const mxArray *mx, const char *field)
{
  insist(mx);
  mxArray *mxPtr = mxGetField(mx, 0, field);
  if(!mxPtr) return 0;
  
  insist(mxIsChar(mxPtr));
  
  char *string = mxArrayToString(mxPtr);
  return string;
}

#endif /* __cplusplus */

#endif /* MATLABUTILS_H */
