#ifndef DEF_H__
#define DEF_H__
#include <float.h>

/* define a real type */
#ifdef HAVE_REAL
  #ifndef HAVEREAL
  #define HAVEREAL HAVE_REAL
  #endif
#endif

#ifndef HAVEREAL
  #define HAVEREAL 1
  typedef double real;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif
