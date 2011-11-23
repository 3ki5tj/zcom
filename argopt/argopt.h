#ifndef INLINE
#define INLINE __inline static
#endif
#include "def.h"
#ifndef ARGOPT_H__
#define ARGOPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  int isopt; /* is option or argument */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  
  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  unsigned flags;
} opt_t;

typedef struct {
  int narg, nopt;
  opt_t *args;
  opt_t *opts;
  const char *prog;
  const char *desc;
  const char *author;
  const struct tm *tm; /* compilation time */
  int version;
  unsigned flags;
  int dum_[4]; /* space holder */
} argopt_t;

#define ARGOPT_MUST     0x0001  /* a mandatory argument or option */
#define ARGOPT_SWITCH   0x0002  /* an option is a switch, fmt "%b" */
#define ARGOPT_SET      0x0004  /* an argument/option is set, fmt starts with "!" */
#define ARGOPT_LONGOPT  0x0010  /* always assume long format, e.g., -maxh */

argopt_t *argopt_open(unsigned flags);
void argopt_close(argopt_t *ao);
#define argopt_regarg(ao, fmt, ptr, desc) argopt_add(ao, NULL, fmt, ptr, desc)
#define argopt_regopt argopt_add
int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc);
void argopt_parse(argopt_t *ao, int argc, char **argv); 

#define argopt_reghelp(ao, sflag) argopt_regopt(ao, sflag, "%b", ao->dum_, "$HELP")
#define argopt_regversion(ao, sflag) argopt_regopt(ao, sflag, "%b", ao->dum_, "$VERSION")

#define argopt_getopt(ao, p) argopt_searchls_(ao->opts, ao->nopt, p)
#define argopt_getarg(ao, p) argopt_searchls_(ao->args, ao->narg, p)
INLINE opt_t *argopt_searchls_(opt_t *ls, int n, const void *p)
 { int i; for (i = 0; i < n; i++) if (ls[i].ptr == p) return ls+i; return NULL; }

#define argopt_set(ao, var) argopt_set_(ao, &var, #var)
INLINE int argopt_set_(argopt_t *ao, const void *p, const char *var)
{ 
   opt_t *a = argopt_getarg(ao, p), *o = argopt_getopt(ao, p);
   die_if(!a && !o, "cannot locate var %s, %p\n", var, p);
   return a ? (a->flags & ARGOPT_SET ? 1 : 0) : (o->flags & ARGOPT_SET ? 1 : 0);
}

#endif

