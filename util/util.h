#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef RESTRICT
#define RESTRICT __restrict
#endif
#ifdef __INTEL_COMPILER
/* operands evaluated in unspecified order */
#pragma warning(disable:981)
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



/* define int16_t/int32_t/int64_t, etc. */
#if (  (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) \
     || defined(__GNUC__) || defined(__INTEL_COMPILER) )
  /* C99 compatible compilers support int64_t etc.
   * but GCC and other compilers has the header even in C89/C90 mode
   * So we need to include more compilers here, see the list on
   * http://sourceforge.net/p/predef/wiki/Compilers/ */
  #include <inttypes.h>
#elif (defined(_MSC_VER) \
      || (defined(__BORLANDC__) && (__BORLANDC__ >= 0x520)))
  /* tested for Visual C++ 6.0 and Borland C++ 5.5 */
  typedef __int8              int8_t;
  typedef __int16             int16_t;
  typedef __int32             int32_t;
  typedef __int64             int64_t;
  typedef unsigned __int8     uint8_t;
  typedef unsigned __int16    uint16_t;
  typedef unsigned __int32    uint32_t;
  typedef unsigned __int64    uint64_t;
#elif defined(__unix__)
  /* a unix compiler is likely to have inttypes.h  */
  #include <inttypes.h>
#else
  /* note the following is a guess, long long is not supported
   * until a later version of visual C++ */
  typedef char                int8_t;
  typedef short               int16_t;
  typedef int                 int32_t;
  typedef long long           int64_t;
  typedef unsigned char       uint8_t;
  typedef unsigned short      uint16_t;
  typedef unsigned            uint32_t;
  typedef unsigned long long  uint64_t;
#endif



/* constant 64-bit integer */
#if defined(_MSC_VER) || defined(__BORLANDC__)
  #define CI64(x) (x ## i64)
  #define CU64(x) (x ## ui64)
#else
  #define CI64(x) (x ## ll)
  #define CU64(x) (x ## ull)
#endif



/* printf() format strings for integers
 * the macros PRId32, PRIu64, etc are defined by C99
 * we write the macros below just in case they are not defined */
#if defined(_MSC_VER) || defined(__BORLANDC__)
  #ifndef PRId32
  #define PRId32 "I32d"
  #endif
  #ifndef SCNd32
  #define SCNd32 "I32d"
  #endif
  #ifndef PRIi32
  #define PRIi32 "I32i"
  #endif
  #ifndef SCNi32
  #define SCNi32 "I32i"
  #endif
  #ifndef PRIu32
  #define PRIu32 "I32u"
  #endif
  #ifndef SCNu32
  #define SCNu32 "I32u"
  #endif
  #ifndef PRIo32
  #define PRIo32 "I32o"
  #endif
  #ifndef SCNo32
  #define SCNo32 "I32o"
  #endif
  #ifndef PRIx32
  #define PRIx32 "I32x"
  #endif
  #ifndef SCNx32
  #define SCNx32 "I32x"
  #endif
  #ifndef PRId64
  #define PRId64 "I64d"
  #endif
  #ifndef SCNd64
  #define SCNd64 "I64d"
  #endif
  #ifndef PRIi64
  #define PRIi64 "I64i"
  #endif
  #ifndef SCNi64
  #define SCNi64 "I64i"
  #endif
  #ifndef PRIu64
  #define PRIu64 "I64u"
  #endif
  #ifndef SCNu64
  #define SCNu64 "I64u"
  #endif
  #ifndef PRIo64
  #define PRIo64 "I64o"
  #endif
  #ifndef SCNo64
  #define SCNo64 "I64o"
  #endif
  #ifndef PRIx64
  #define PRIx64 "I64x"
  #endif
  #ifndef SCNx64
  #define SCNx64 "I64x"
  #endif
#else
  #ifndef PRId32
  #define PRId32 "d"
  #endif
  #ifndef SCNd32
  #define SCNd32 "d"
  #endif
  #ifndef PRIi32
  #define PRIi32 "i"
  #endif
  #ifndef SCNi32
  #define SCNi32 "i"
  #endif
  #ifndef PRIu32
  #define PRIu32 "u"
  #endif
  #ifndef SCNu32
  #define SCNu32 "u"
  #endif
  #ifndef PRIo32
  #define PRIo32 "o"
  #endif
  #ifndef SCNo32
  #define SCNo32 "o"
  #endif
  #ifndef PRIx32
  #define PRIx32 "x"
  #endif
  #ifndef SCNx32
  #define SCNx32 "x"
  #endif
  #ifndef PRId64
  #define PRId64 "lld"
  #endif
  #ifndef SCNd64
  #define SCNd64 "lld"
  #endif
  #ifndef PRIi64
  #define PRIi64 "lli"
  #endif
  #ifndef SCNi64
  #define SCNi64 "lli"
  #endif
  #ifndef PRIu64
  #define PRIu64 "llu"
  #endif
  #ifndef SCNu64
  #define SCNu64 "llu"
  #endif
  #ifndef PRIo64
  #define PRIo64 "llo"
  #endif
  #ifndef SCNo64
  #define SCNo64 "llo"
  #endif
  #ifndef PRIx64
  #define PRIx64 "llx"
  #endif
  #ifndef SCNx64
  #define SCNx64 "llx"
  #endif
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



#ifndef xnew
#define xnew(x, n) { \
  size_t num_ = (size_t) (n); \
  die_if (num_ <= 0, \
    "cannot allocate %d objects for %s\n", (int) num_, #x); \
  die_if ((x = calloc(num_, sizeof(*(x)))) == NULL, \
    "no memory for %s x %d\n", #x, (int) num_); }
#endif



#ifndef xrenew
#define xrenew(x, n) { \
  size_t num_ = (size_t) (n); \
  die_if (num_ <= 0, \
    "cannot allocate %d objects for %s\n", (int) num_, #x); \
  die_if ((x = realloc(x, (num_)*sizeof(*(x)))) == NULL, \
    "no memory for %s x %d\n", #x, (int) num_); }
#endif



#define xfopen(fp, fn, fmt, err) \
  if ((fp = fopen(fn, fmt)) == NULL) { \
    fprintf(stderr, "cannot open file %s with mode [%s]\n", fn, fmt); \
    err; }



/* check if file `fn' exists */
INLINE int fexists(const char *fn)
{
  FILE *fp;
  if ((fp = fopen(fn, "r")) == NULL) return 0;
  else { fclose(fp); return 1; }
}



/* copy file */
INLINE int copyfile(const char *fninp, const char *fnout)
{
  FILE *fpinp, *fpout;
#ifndef COPYFILE_BUFSZ
#define COPYFILE_BUFSZ (64*1024)
#endif
  unsigned char buf[COPYFILE_BUFSZ];
  size_t sz, tot = 0;

  if ((fpinp = fopen(fninp, "rb")) == NULL) {
    fprintf(stderr, "copyfile: cannot read file %s\n", fninp);
    return -1;
  }
  if ((fpout = fopen(fnout, "wb")) == NULL) {
    fprintf(stderr, "copyfile: cannot write file %s\n", fnout);
    fclose(fpout);
    return -2;
  }
  while ((sz = fread(buf, sizeof(buf[1]), COPYFILE_BUFSZ, fpinp)) != 0) {
    tot += sz;
    /* note: sz may differ from COPYFILE_BUFSZ */
    if (sz != fwrite(buf, sizeof(buf[1]), sz, fpout))
      fprintf(stderr, "copyfile: error writing %s, byte %.0f\n", fnout, 1.*tot);
    if ( feof(fpinp) ) break;
  }
  fclose(fpinp);
  fclose(fpout);
  return 0;
}



/* swap two variables */
#define xtpswap(tp, x, y) { tp dum_; dum_ = (x); (x) = (y); (y) = dum_; }

#define intswap(x, y) xtpswap(int, x, y)

#define dblswap(x, y) xtpswap(double, x, y)

#define realswap(x, y) xtpswap(real, x, y)



INLINE int intmax(int x, int y) { return x > y ? x : y; }
INLINE int intmin(int x, int y) { return x < y ? x : y; }



/* confine x within [xmin, xmax] */
INLINE int intconfine(int x, int xmin, int xmax)
  { return x < xmin ? xmin : x > xmax ? xmax : x; }



INLINE int intsqr(int x) { return x * x; }



/* get the pair index from 0 to n*(n - 1)/2 - 1 */
INLINE int getpairindex(int i, int j, int n)
{
  die_if (i < 0 || i >= n || j < 0 || j >= n || i == j,
      "bad index error i %d, j %d, n %d\n", i, j, n);
  if (i > j) { int i1 = i; i = j; j = i1; }
  return n*i - (i + 1)*(i + 2)/2 + j;
}

/* return individual indices for a given pair index */
INLINE int parsepairindex(int id, int n, int *i, int *j)
{
  int i1, n1;
  die_if (id < 0 || id >= n*(n - 1)/2, "index %d too large for n %d\n", id, n);
  for (i1 = n - 2; i1 >= 0; i1--) {
    if (id >= (n1 = i1*n - i1*(i1 + 1)/2)) {
      *i = i1;
      *j = id - n1 + i1 + 1;
      return 0;
    }
  }
  return 1;
}



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
  return x*sqrt(1 + t*t);
}



/* round x to a multiple dx  */
INLINE double dblround(double x, double dx)
{
  if (x*dx > 0) return dx * (int)(x/dx + (.5 - DBL_EPSILON));
  else return -dx * (int)(-x/dx + (.5 - DBL_EPSILON));
}



/* convert to double to integer */
INLINE int dbl2int(double x)
{
  return (int) ((x < 0) ? (x - .5) : (x + .5));
}



INLINE void dblcleararr(double *x, int n)
{
  int i; for (i = 0; i < n; i++) x[i] = 0.0;
}



#ifdef HAVEREAL
INLINE real realmax(real x, real y) { return x > y ? x : y; }
INLINE real realmin(real x, real y) { return x < y ? x : y; }
/* confine x within [xmin, xmax] */
INLINE real realconfine(real x, real xmin, real xmax)
{ return x < xmin ? xmin : x > xmax ? xmax : x; }
#endif



#ifndef LNADD_DEFINED
#define LNADD_DEFINED
#define LN_BIG 50.0

/* log(exp(a) + exp(b)) */
INLINE double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > LN_BIG) ? a : a + log(1 + exp(-c));
}

/* log(exp(a) - exp(b)), only works for a > b */
INLINE double lndif(double a, double b)
{
  double c;
  die_if (a < b, "lndif: %g < %g\n", a, b);
  return ((c = a - b) > LN_BIG) ? a : a + log(1 - exp(-c));
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
    else memmove(s, p, strlen(p) + 1);
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
 * len non-blank characters, i.e., s[len] = '\0' for the longest output */
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
  return cs - ct;
}



#endif

