
/* $Id$ */

#include <unistd.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include "die.h"

void die_at(const char *s, const char *file, int line)
{
  MatlabCrashLoc(s, file, line);
}

void die(const char *s)
{
  cout << s << "\n" << endl;
}



