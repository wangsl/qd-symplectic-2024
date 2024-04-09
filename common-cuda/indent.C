

/* $Id$ */

#include "indent.h"

static char spaces[21] = "                    ";
static int spi = 20;

void IndentPush()
{
  spi -= 2;
}

void IndentPop()
{
  spi += 2;
}

const char *Indent()
{
  int i = spi;
  if (i > 20)
    i = 20;
  else if (i < 0)
    i = 0;
  return spaces + i;
}

