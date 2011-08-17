#define ZCINLINE __inline static

#ifndef UTIL_H__
#define UTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#ifndef xnew
#define xnew(x, n) \
  if ((n) <= 0) { \
    fprintf(stderr, "cannot allocate %d objects for %s\n", (int) (n), #x); \
    exit(1); \
  } else if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); } 
#endif

#ifndef xrenew
#define xrenew(x, n) \
  if ((n) <= 0) { \
    fprintf(stderr, "cannot allocate %d objects for %s\n", (int) (n), #x); \
    exit(1); \
  } else if ((x = realloc(x, (n)*sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); } 
#endif

/* error handling */
void perrmsg__(const char *file, int line, const char *why,
    int err, const char *fmt, va_list args);

#ifdef ZCHAVEVAM

ZCINLINE void perrmsg_(const char *file, int line, const char *why,
    int cond, int err, const char *fmt, ...)
{
  va_list args;
 
  if (cond) {
    va_start(args, fmt);
    perrmsg__(file, line, why, err, fmt, args);
    va_end(args);
    if (err) exit(1);
  }
}
#define die_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, 1, fmt, ## __VA_ARGS__)
#define msg_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, 0, fmt, ## __VA_ARGS__)
#define fatal(fmt, ...)  die_if(1, fmt, ## __VA_ARGS__)

#else /* ZCHAVEVAM */

#define PERRMSG__(c, x) {                     \
  va_list args;                               \
  if ((#c[0] == '1' && #c[1] == '\0') || c) { \
    va_start(args, fmt);                      \
    perrmsg__(NULL, -1, NULL, x, fmt, args);  \
    va_end(args);                             \
    if (#x[0] == '1') exit(1);                \
  } }
ZCINLINE void die_if(int cond, const char *fmt, ...) PERRMSG__(cond, 1)
#ifdef USE_MSG_IF
ZCINLINE void msg_if(int cond, const char *fmt, ...) PERRMSG__(cond, 0)
#endif
#ifdef USE_FATAL
void fatal(const char *fmt, ...) PERRMSG__(1, 1)
#endif
#undef PERRMSG__

#endif /* ZCHAVEVAM */

ZCINLINE int fexists(const char *fn)
{
  FILE *fp; 
  if ((fp = fopen(fn, "r")) == NULL) return 0;
  else { fclose(fp); return 1; }
} 

/* sqrt(x*x + y*y) */
ZCINLINE double hypotn(double x, double y)
{
  double t;
  x = fabs(x);
  y = fabs(y);
  if (x <= 0.) return y;
  else if (y <= 0.) return x;
  if (x < y) t = x, x = y, y = t;
  t = y/x;
  return x*sqrt(1+t*t);
}

/* round x to a multiple dx  */
ZCINLINE double dblround(double x, double dx)
{
  if (x > 0) return dx * (int)(x/dx+.5-1e-14);
  else return -dx * (int)(-x/dx+.5-1e-14);
}

ZCINLINE double dblsqr(double x) { return x*x; }

#ifndef LNADD_DEFINED
#define LNADD_DEFINED
#define LN_BIG 50.0

/* log(exp(a) + exp(b)) */
ZCINLINE double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a-b) > LN_BIG) ? a : a + log(1 + exp(-c));
}

/* log(exp(a)-exp(b)), only works for a>b */
ZCINLINE double lnmin(double a, double b)
{
  double c;
  die_if (a < b, "lnmin: %g < %g\n", a, b);
  return ((c = a-b) > LN_BIG) ? a : a + log(1 - exp(-c));
}

/* log(exp(a)+b) */
ZCINLINE double lnaddn(double a, double b)
{
  return (a > LN_BIG) ? a : a + log(1 + b*exp(-a));
}

#undef LN_BIG
#endif

#endif

