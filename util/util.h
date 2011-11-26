#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef UTIL_H__
#define UTIL_H__
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#ifndef xnew
#define xnew(x, n) \
  if (#n[0] != '1' && (n) <= 0) { \
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

/* print an error message */
INLINE void perrmsg__(const char *file, int line, const char *why,
    const char *fmt, va_list args)
{
  fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  if (fmt[strlen(fmt) - 1] != '\n')
    fprintf(stderr, "\n"); /* add a new line if needed */
  if (file != NULL) fprintf(stderr, "file: %s\n", file);
  if (line >= 0) fprintf(stderr, "line: %d\n", line);
  if (why != NULL && strcmp(why, "1") != 0) 
    fprintf(stderr, "cond: %s\n", why);
}

#ifdef HAVEVAM

INLINE void perrmsg_(const char *file, int line, const char *why,
    int cond, const char *fmt, ...)
{
  if (cond) {
    va_list args;
    va_start(args, fmt);
    perrmsg__(file, line, why, fmt, args);
    va_end(args);
    exit(1);
  }
}

#define die_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, fmt, ## __VA_ARGS__)
#define fatal(fmt, ...)  die_if(1, fmt, ## __VA_ARGS__)

#else /* !HAVEVAM */

#define PERRMSG__(c) {                        \
  if ((#c[0] == '1' && #c[1] == '\0') || c) { \
    va_list args;                             \
    va_start(args, fmt);                      \
    perrmsg__(NULL, -1, NULL, fmt, args);     \
    va_end(args);                             \
    exit(1);                                  \
  } }
INLINE void die_if(int cond, const char *fmt, ...) PERRMSG__(cond)
INLINE void fatal(const char *fmt, ...) PERRMSG__(1)
#undef PERRMSG__

#endif /* HAVEVAM */

#define xfopen(fp, fn, fmt, err) \
  if ((fp = fopen(fn, fmt)) == NULL) { \
    fprintf(stderr, "cannot open file %s\n", fn); err; }

INLINE int fexists(const char *fn)
{
  FILE *fp; 
  if ((fp = fopen(fn, "r")) == NULL) return 0;
  else { fclose(fp); return 1; }
} 

/* swap two variables */
#ifndef xtpswap
#define xtpswap(tp, x, y) { tp dum_; dum_ = (x); (x) = (y); (y) = dum_; }
#endif

#ifndef intswap
#define intswap(x, y) xtpswap(int, x, y)
#endif

#ifndef dblswap
#define dblswap(x, y) xtpswap(double, x, y)
#endif

INLINE int intmax(int x, int y) { return x > y ? x : y; }
INLINE int intmin(int x, int y) { return x < y ? x : y; }
/* confine x within [xmin, xmax] */
INLINE int intconfine(int x, int xmin, int xmax)
  { return x < xmin ? xmin : x > xmax ? xmax : x; }

INLINE int intsqr(int x) { return x * x; }

INLINE double dblmax(double x, double y) { return x > y ? x : y; }
INLINE double dblmin(double x, double y) { return x < y ? x : y; }
/* confine x within [xmin, xmax] */
INLINE double dblconfine(double x, double xmin, double xmax)
  { return x < xmin ? xmin : x > xmax ? xmax : x; }

INLINE double dblsqr(double x) { return x * x; }

/* sqrt(x*x + y*y) */
INLINE double dblhypot(double x, double y)
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
INLINE double dblround(double x, double dx)
{
  if (x*dx > 0) return dx * (int)(x/dx + (.5 - DBL_EPSILON));
  else return -dx * (int)(-x/dx + (.5 - DBL_EPSILON));
}

INLINE void dblcleararr(double *x, int n)
  { int i; for (i = 0; i < n; i++) x[i] = 0.0; }

#ifndef LNADD_DEFINED
#define LNADD_DEFINED
#define LN_BIG 50.0

/* log(exp(a) + exp(b)) */
INLINE double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a-b) > LN_BIG) ? a : a + log(1 + exp(-c));
}

/* log(exp(a) - exp(b)), only works for a > b */
INLINE double lndif(double a, double b)
{
  double c;
  die_if (a < b, "lndif: %g < %g\n", a, b);
  return ((c = a-b) > LN_BIG) ? a : a + log(1 - exp(-c));
}

/* log(exp(a)+b) */
INLINE double lnaddn(double a, double b)
{
  return (a > LN_BIG) ? a : a + log(1 + b*exp(-a));
}

#undef LN_BIG
#endif  /* LNADD_DEFINED */

#define cisalnum(c)   isalnum((unsigned char)(c))
#define cisalpha(c)   isalpha((unsigned char)(c))
#define cisdigit(c)   isdigit((unsigned char)(c))
#define cisxdigit(c)  isxdigit((unsigned char)(c))
#define cisprint(c)   isprint((unsigned char)(c))
#define cisspace(c)   isspace((unsigned char)(c))
#define cislower(c)   islower((unsigned char)(c))
#define cisupper(c)   isupper((unsigned char)(c))
#define ctolower(c)   (char) tolower((unsigned char)(c))
#define ctoupper(c)   (char) toupper((unsigned char)(c))

/* string manipulation */
#define ZSTR_XSPACEL  0x0001
#define ZSTR_XSPACER  0x0002
#define ZSTR_XSPACE   (ZSTR_XSPACEL|ZSTR_XSPACER)
#define ZSTR_COPY     0x0004
#define ZSTR_CAT      0x0008
#define ZSTR_CASE     0x0100
#define ZSTR_UPPER_   0x0200
#define ZSTR_UPPER    (ZSTR_CASE|ZSTR_UPPER_)
#define ZSTR_LOWER    ZSTR_CASE

/* remove leading and trailing spaces */
#define strip(s)  stripx(s, ZSTR_XSPACE)
#define lstrip(s) stripx(s, ZSTR_XSPACEL)
#define rstrip(s) stripx(s, ZSTR_XSPACER)
INLINE char *stripx(char *s, unsigned flags)
{
  char *p;

  if (flags & ZSTR_XSPACEL) { /* remove leading spaces */
    for (p = s; cisspace(*p); p++) ;
    if (*p == '\0') *s = '\0';
    else memmove(s, p, strlen(p)+1);
  }
  if (flags & ZSTR_XSPACER) /* remove trailing spaces */
    for (p = s + strlen(s) - 1; p >= s && cisspace(*p); p--)
      *p = '\0';
  return s;
}

/* in the follows, size_s means the buffer size of s, i.e., sizeof(s) for static strings */
/* copy the string and convert it to upper/lower case */
#define strcpy2u(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY|ZSTR_UPPER)
#define strcpy2l(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY|ZSTR_LOWER)
#define strcpy_sf(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY)
#define substr(s, t, start, len) strcnv(s, t+start, len, ZSTR_COPY)
/* concatenate strings, the last parameter is the buffer size of s,
 * unlike strncat(), in which it's the number of characters from *t* to be copied.  */
#define strcat_sf(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_CAT)
/* safely copy/cat strings with case conversion
 * unlike strncpy(), s is always null-terminated on return: it copies at most
 * len nonblank characters, i.e., s[len] = '\0' for the longest output */
INLINE char *strcnv(char *s, const char *t, size_t len, unsigned flags)
{
  size_t i = 0, j;
  unsigned docase = flags & ZSTR_CASE, up = flags & ZSTR_UPPER_;

  if (len == 0 || s == NULL || t == NULL) return s;
  if (flags & ZSTR_CAT) while(s[i]) i++;
  for (j = 0; i < len; i++, j++) {
    if (docase && t[j]) {
      if (up) s[i] = (char) (unsigned char) toupper((unsigned char) t[j]);
      else    s[i] = (char) (unsigned char) tolower((unsigned char) t[j]);
    } else s[i] = t[j];
    if (t[j] == 0) break;
  }
  if (i == len) s[i] = '\0';
  if (flags & ZSTR_XSPACE) stripx(s, flags); /* call strip */
  return s;
}

/* compare strings without case */
#define strcmpnc(s, t) strncmpnc(s, t, -1)
INLINE int strncmpnc(const char *s, const char *t, int n)
{
  int i, cs, ct;

  if (s == NULL || t == NULL) return 0;
  for (i = 0; ; i++) {
    if (i >= n) return 0;
    cs = s[i];
    ct = t[i];
    if (cs == 0 || ct == 0) break;
    cs = toupper( (unsigned char) cs );
    ct = toupper( (unsigned char) ct );
    if (cs != ct) break;
  }
  return cs-ct;
}

#endif

