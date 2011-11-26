#include "def.h"
#include "opt.h"
#ifndef ARGOPT_H__
#define ARGOPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

#define ARGOPT_MUST     OPT_MUST    /* mandatory argument or option, format starts with ! */
#define ARGOPT_SWITCH   OPT_SWITCH  /* format "%b" */
#define ARGOPT_SET      OPT_SET
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

 /* test if argument/option is set */
#define argopt_set(ao, var) argopt_set_(ao, &var, #var)
INLINE int argopt_set_(argopt_t *ao, const void *p, const char *var)
{ 
   opt_t *a = argopt_getarg(ao, p), *o = argopt_getopt(ao, p);
   die_if(!a && !o, "cannot locate var %s, %p\n", var, p);
   return a ? (a->flags & OPT_SET ? 1 : 0) : (o->flags & OPT_SET ? 1 : 0);
}

#endif

