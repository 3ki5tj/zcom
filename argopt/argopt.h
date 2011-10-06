#ifndef ZCINLINE
#define ZCINLINE __inline static
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
  int dum_;
} argopt_t;

#define ARGOPT_MUST     0x0001  /* a mandatory argument or option */
#define ARGOPT_SWITCH   0x0002  /* an option is a switch, fmt "%b" */
#define ARGOPT_SET      0x0004  /* an argument/option is set, fmt starts with "!" */
#define ARGOPT_LONGOPT  0x0010  /* always assume long format, e.g., -maxh */

argopt_t *argopt_open(unsigned flags);
void argopt_close(argopt_t *ao);
int argopt_regarg(argopt_t *ao, const char *fmt, void *ptr, const char *desc);
int argopt_regopt(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc);
void argopt_parse(argopt_t *ao, int argc, char **argv); 

#define argopt_regopt_help(ao, sflag) argopt_regopt(ao, sflag, NULL, &ao->dum_, "$HELP")
#define argopt_regopt_version(ao, sflag) argopt_regopt(ao, sflag, NULL, &ao->dum_, "$VERSION")

ZCINLINE opt_t *argopt_getopt(argopt_t *ao, const void *p)
 { int i; for (i = 0; i < ao->nopt; i++) if (ao->opts[i].ptr == p) return ao->opts+i; return NULL; }
ZCINLINE opt_t *argopt_getarg(argopt_t *ao, const void *p)
 { int i; for (i = 0; i < ao->narg; i++) if (ao->args[i].ptr == p) return ao->args+i; return NULL; }

#define argopt_set(ao, var) argopt_set_(ao, &var, #var)
ZCINLINE int argopt_set_(argopt_t *ao, const void *p, const char *var)
{ 
   opt_t *a = argopt_getarg(ao, p), *o = argopt_getopt(ao, p);
   die_if(!a && !o, "cannot locate var %s, %p\n", var, p);
   return a ? (a->flags & ARGOPT_SET ? 1 : 0) : (o->flags & ARGOPT_SET ? 1 : 0);
}

#endif

