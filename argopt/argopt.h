#ifndef ARGOPT_H__
#define ARGOPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  const char *val; /* raw string from the command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  void *ptr; /* address to the target variable */
} arg_t;

typedef struct {
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *val; /* raw string from command line */
  const char *desc; /* description */
  int hasarg;  /* needs an additional argument */
  int set; /* if the option has been set or not */
  const char *fmt; /* sscanf format */
  void *ptr; /* address to the target variable */
} opt_t;

typedef struct {
  int narg, nopt;
  arg_t *args;
  opt_t *opts;
  const char *prog;
} argopt_t;

#endif
