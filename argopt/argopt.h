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
  int must; /* must be given */
  int hasarg;  /* needs an additional argument */
  int set; /* if the option has been set or not */
} opt_t;

typedef struct {
  int narg, nopt;
  opt_t *args;
  opt_t *opts;
  const char *prog;
  const char *desc;
  const char *author;
  int version;
  const struct tm *tm; /* compilation time */
  int dum_;
} argopt_t;

argopt_t *argopt_open(int version, const char *desc, const char *author);
void argopt_close(argopt_t *ao);
opt_t *argopt_regarg(argopt_t *ao, int must, 
    const char *fmt, void *ptr, const char *desc);
opt_t *argopt_regopt(argopt_t *ao, const char *sflag, int hasarg,
    const char *fmt, void *ptr, const char *desc);
void argopt_parse(argopt_t *ao, int argc, const char **argv); 

#define argopt_regopt_help(ao, sflag) argopt_regopt(ao, sflag, 0, NULL, &ao->dum_, "$HELP")
#define argopt_regopt_version(ao, sflag) argopt_regopt(ao, sflag, 0, NULL, &ao->dum_, "$VERSION")

#endif

