#ifndef DEF_H__
#define DEF_H__

/* define a real type */
#ifdef HAVE_REAL
  #ifndef ZCHAVEREAL
  #define ZCHAVEREAL HAVE_REAL
  #endif
#endif

#ifndef ZCHAVEREAL
  #define ZCHAVEREAL 1
  typedef double real;
#endif

#endif
