#include "def.h"
#include "opt.h"
#ifndef ARGOPT_H__
#define ARGOPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  int nopt;
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
#define argopt_reghelp argopt_addhelp
#define argopt_regversion argopt_addversion
int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc);
void argopt_parse(argopt_t *ao, int argc, char **argv);
INLINE void argopt_dump(const argopt_t *ao);

#define argopt_addhelp(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$HELP")
#define argopt_addversion(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$VERSION")

#define argopt_getopt(ao, p) opt_find(ao->opts, ao->nopt, p)
#define argopt_getarg argopt_getopt

/* test if argument/option is explicitly set */
#define argopt_set(ao, var) opt_isset(ao->opts, ao->nopt, &var, #var)

#endif

