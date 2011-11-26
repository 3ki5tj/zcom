#ifndef INLINE
#define INLINE __inline static
#endif
#include "def.h"
#include "ss.c"
#ifndef OPT_H__
#define OPT_H__
#include <stdio.h>

/* option either from arguments or configuration */
typedef struct {
  int isopt; /* is option or argument */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *key; /* key */
  
  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  unsigned flags;
} opt_t;

#define OPT_MUST     0x0001  /* a mandatory argument or option */
#define OPT_SWITCH   0x0002  /* an option is a switch */
#define OPT_SET      0x0004  /* an argument/option is set */

/* translate string values to actual ones through sscanf() */
INLINE int opt_getval(opt_t *o)
{
  const char *fmt = o->fmt;
  
  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
    *(const char **)o->ptr = o->val;
  } else if (strcmp(fmt, "%s") == 0) {
    sscpy( *(char **)o->ptr, o->val);
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

/* set properties of an option: fmt = "%b" for a switch */
INLINE void opt_set(opt_t *o, const char *sflag, const char *key,
    const char *fmt, void *ptr, const char *desc)
{
  o->ch = '\0';
  if (key) {
    o->isopt = 2;
  } else if (sflag) { /* option */
    o->isopt = 1;
    o->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    o->isopt = 0;
  }
  o->sflag = sflag;
  o->key = key;
  o->flags = 0;
  die_if (ptr == NULL, "null pass to argopt with %s: %s\n", sflag, desc);
  o->ptr = ptr;
  if (fmt == NULL) fmt = "";
  if (fmt[0] == '!') {
    fmt++;
    o->flags |= OPT_MUST;
  }
  if (strcmp(fmt, "%b") == 0) {
    fmt = "%d";
    o->flags |= OPT_SWITCH;
  }
  o->fmt = fmt;
  o->pfmt = NULL;
  o->desc = desc;
}

/* print the value of o->ptr */
INLINE void opt_printptr(opt_t *o)
{
  const char *fmt;

  for (fmt = o->fmt; *fmt && *fmt != '%'; fmt++) ;
#define ELIF_PF_(fm, fmp, type) else if (strcmp(fmt, fm) == 0) printf((o->pfmt ? o->pfmt : fmp), *(type *)o->ptr)
  if (fmt == NULL || *fmt == '\0') printf("%s", (*(char **)o->ptr) ? (*(char **)o->ptr) : "NULL");
  ELIF_PF_("%b", "%d", int); /* switch */
  ELIF_PF_("%d", "%d", int);
  ELIF_PF_("%u", "%u", unsigned);
  ELIF_PF_("%x", "0x%x", unsigned);
  ELIF_PF_("%ld", "%ld", long);
  ELIF_PF_("%lu", "%lu", unsigned long);
  ELIF_PF_("%lx", "0x%lx", unsigned long);
#if 0  /* C99 only */
  ELIF_PF_("%lld", "%lld", long long);
  ELIF_PF_("%llu", "%llu", unsigned long long);
  ELIF_PF_("%llx", "0x%llx", unsigned long long);
#endif
  ELIF_PF_("%f", "%g", float);
  ELIF_PF_("%lf", "%g", double);
  ELIF_PF_("%r", "%g", real);
  else printf("unknown %s-->%%d: %d\n", fmt, *(int *)o->ptr);
#undef ELIF_PF_  
}

#endif
