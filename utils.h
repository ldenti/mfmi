#ifndef UTILS_H_
#define UTILS_H_

#include <sys/time.h>

static inline double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}

#endif
