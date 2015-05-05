#include "util.h"
#include "ss.h"
#ifndef OPT_H__
#define OPT_H__



enum { OPT_ARGUMENT, OPT_OPTION, OPT_CFG, OPT_COUNT };



/* option either from arguments or configuration */
typedef struct {
  int isopt; /* one of OPT_xxx values */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *key; /* key, for cfg files as in `key = val' */

  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  int initval; /* initial value, for a switch option */
  const char **sarr; /* array of string values, for a list option */
  int scnt; /* length of the array, for a list option */
  unsigned flags;
} opt_t;



/* support __float128 for GCC
 * Intel compiler defines __GNUC__, but it does not have __float128 */
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) \
  && ( (__GNUC__ == 4 && __GNUC_MINOR__ >= 6) || __GNUC__ > 4 )
  #define HAVEFLOAT128 1
  #include <quadmath.h>
  /* ignore warnings for "%Qf"
   * assume `diagnostic push' is available in this case */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wformat"
  #pragma GCC diagnostic ignored "-Wformat-extra-args"
  #pragma GCC diagnostic ignored "-Wlong-long"
#else
  #define HAVEFLOAT128 0
#endif



/* return the index of string from a predefined array
 * using fuzzy string comparison */
INLINE int opt_select(const char *s, const char **sarr, int n)
{
  int i;
  char ls[1024];

  for ( i = 0; i < n; i++ )
    if ( strcmpfuzzy(sarr[i], s) == 0 )
      return i;
  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i >= 0 && i < n ) return i;
  }
  strjoin(ls, sizeof ls, sarr, n, ", ");
  fprintf(stderr, "Error: cannot select %s from the array of %d items: %s\n",
      s, n, ls);
  exit(1);
  return 0;
}



#define OPT_MUST     0x0001  /* a mandatory argument or option */
#define OPT_SWITCH   0x0002  /* an option is a switch */
#define OPT_SET      0x0004  /* an argument/option is set */

/* translate string value in `o->val' into
 * actual ones through sscanf(), etc */
INLINE int opt_getval(opt_t *o)
{
  const char *fmt = o->fmt;

  /* for a string argument, it can be obtained from NULL or "%s"
   * NULL: string memory from command-line
   * "%s":  string memory from sscpy()
   * if the string is read-only, like a output file name,
   * but the string is to be changed, use "%s" */
  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
    *((const char **) o->ptr) = o->val;
  } else if (strcmp(fmt, "%s") == 0) { /* copy the string */
    sscpy( *((char **) o->ptr), o->val);
  } else if (strcmpnc(fmt, "%list") == 0
          || strcmpnc(fmt, "%enum") == 0) {
    *((int *) o->ptr) = opt_select(o->val, o->sarr, o->scnt);
  } else if (strcmp(fmt, "%b") == 0) { /* switch */
    /* switch the default value */
    if (o->flags & OPT_SET) return !o->initval;
    else return o->initval;
  } else if (strcmp(fmt, "%+") == 0) { /* incremental, like -vv */
    *((int *) o->ptr) += 1;
  } else if (strcmp(fmt, "%-") == 0) { /* decremental, like -vv */
    *((int *) o->ptr) -= 1;
#if HAVEFLOAT128
  } else if (strcmp(fmt, "%Qf") == 0) {
#if defined(QUAD) || defined(F128)
    *((__float128 *) o->ptr) = strtoflt128(o->val, NULL);
#else /* intrinsic hook, requires no -lquadmath if unnecessary */
    sscanf(o->val, "%Qf", (__float128 *) o->ptr);
#endif /* defined(QUAD) || defined(F128) */
#endif
  } else { /* call sscanf */
    if (strcmp(fmt, "%r") == 0) /* real */
      fmt = (sizeof(real) == sizeof(float)) ? "%f" : "%lf";
    if (1 != sscanf(o->val, fmt, o->ptr)) {
      fprintf(stderr, "Error: unable to convert a value for [%s] as fmt [%s], raw string: [%s]\n",
          o->desc, fmt, o->val);
      return 1;
    }
  }
  return 0;
}



/* register an option
 *
 * for a configure entry, set `key' and leave `sflag' = NULL
 * for a command-line option, set `sflag' and leave `key' = NULL
 * `fmt' is the sscanf() format string
 * `*ptr' is the target variable
 * `fmt' can "%b" for a switch (like an command-line option "-v")
 * `fmt' can have a prefix `!' to mean a mandatory option
 * both NULL and "%s" of `fmt' mean string values, the type of
 *  `ptr' should be `char **', the difference is that `*ptr'
 *  is directly assigned to `o->val' during opt_getval() in the
 *  former case, but extra memory is allocated to copy `o->val'
 *  in the latter case */
#define opt_set(o, sflag, key, fmt, ptr, desc) \
  opt_setx(o, sflag, key, fmt, ptr, desc, NULL, 0)

INLINE void opt_setx(opt_t *o, const char *sflag, const char *key,
    const char *fmt, void *ptr, const char *desc,
    const char **sarr, int scnt)
{
  o->ch = '\0';
  if (key) { /* cfg file `key = val', not a command-line argument */
    o->isopt = OPT_CFG;
  } else if (sflag) { /* option */
    o->isopt = OPT_OPTION;
    o->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    o->isopt = OPT_ARGUMENT;
  }
  o->sflag = sflag;
  o->key = key;
  o->flags = 0;
  die_if (ptr == NULL, "null pass to opt with %s: %s\n", sflag, desc);
  o->ptr = ptr;
  if (fmt == NULL) fmt = "";
  if (fmt[0] == '!') {
    fmt++;
    o->flags |= OPT_MUST;
  }
  die_if (fmt[0] != '\0' && fmt[0] != '%',
      "unknown format (missing `%%') flag `%s\', fmt `%s', description: %s\n",
      sflag, fmt, desc);
  if (strcmp(fmt, "%b") == 0) {
    o->flags |= OPT_SWITCH;
    o->initval = *((int *) ptr); /* save the initial value */
  }
  o->fmt = fmt;
  o->pfmt = NULL;
  o->desc = desc;
  o->sarr = sarr;
  o->scnt = scnt;
}



/* print the value of o->ptr */
#define opt_printptr(o) opt_fprintptr(stderr, o)
INLINE void opt_fprintptr(FILE *fp, opt_t *o)
{
  const char *fmt;

  for (fmt = o->fmt; *fmt && *fmt != '%'; fmt++) ;

#define ELIF_PF_(fm, fmp, type) \
  else if (strcmp(fmt, fm) == 0) \
  { fprintf(fp, (o->pfmt ? o->pfmt : fmp), *(type *)o->ptr); }

  if (fmt == NULL || *fmt == '\0' || strcmp(fmt, "%s") == 0) {
    fprintf(fp, "%s", (*(char **) o->ptr) ? (*(char **) o->ptr) : "NULL");
  } else if (strcmp(fmt, "%b") == 0) {
    fprintf(fp, "%s", (*(int *)o->ptr) ? "true" : "false");
  } else if (strcmpfuzzy(fmt, "%list") == 0
          || strcmpfuzzy(fmt, "%enum") == 0) {
    int ival = *((int *) o->ptr);
    fprintf(fp, "%d (%s)", ival, o->sarr[ival]);
  }
  ELIF_PF_("%+", "%d", int)
  ELIF_PF_("%-", "%d", int)
  ELIF_PF_("%d", "%d", int)
  ELIF_PF_("%u", "%u", unsigned)
  ELIF_PF_("%x", "0x%x", unsigned)
  ELIF_PF_("%ld", "%ld", long)
  ELIF_PF_("%lo", "%lo", long)
  ELIF_PF_("%lu", "%lu", unsigned long)
  ELIF_PF_("%lx", "0x%lx", unsigned long)
#if defined(HAVELONGLONG) /* C99 or GCC extension */
  ELIF_PF_("%lld", "%lld", long long)
  ELIF_PF_("%llo", "%llo", long long)
  ELIF_PF_("%llu", "%llu", unsigned long long)
  ELIF_PF_("%llx", "0x%llx", unsigned long long)
#endif
  ELIF_PF_("%f", "%g", float)
  ELIF_PF_("%lf", "%g", double)
  ELIF_PF_("%Lf", "%Lg", long double)
#if HAVEFLOAT128
#if defined(QUAD) || defined(F128)
  else if (strcmp(fmt, "%Qf") == 0) {
    char buf[256];
    quadmath_snprintf(buf, sizeof buf,
        (o->pfmt ? o->pfmt : "%Qg"), *(__float128 *)o->ptr);
    fprintf(fp, "%s", buf);
  }
#else /* intrinsic hook, requires no -lquadmath if unnecessary */
  ELIF_PF_("%Qf", "%Qg", __float128)
#endif /* defined(QUAD) || defined(F128) */
#endif /* HAVEFLOAT128 */
  ELIF_PF_("%r", "%g", real)
  else {
    fprintf(fp, "unknown %s-->%%d: %d", fmt, *(int *) o->ptr);
  }

#undef ELIF_PF_
}



/* search an option list, return an option whose variable address is p */
INLINE opt_t *opt_find(opt_t *ls, int n, const void *p)
{
  int i;

  for (i = 0; i < n; i++) if (ls[i].ptr == p) return ls + i;
  return NULL;
}



/* search an option list to see if an option is explicitly set */
INLINE int opt_isset(opt_t *ls, int n, const void *p, const char *var)
{
  opt_t *o = opt_find(ls, n, p);
  die_if (!o, "cannot find var %s, ptr %p\n", var, p);
  return o->flags & OPT_SET ? 1 : 0;
}



#if HAVEFLOAT128
  #pragma GCC diagnostic pop
#endif



#endif
