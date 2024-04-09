
#ifndef GEN_UTILS_H
#define GEN_UTILS_H

#include <sys/time.h>
#include "cudaUtils.h"

inline double get_time_now_in_secs()
{
  struct timeval t;
  insist(!gettimeofday(&t, NULL));
  return t.tv_sec + t.tv_usec*1.0e-6;
}

#endif /* GEN_UTILS_H */
