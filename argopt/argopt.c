#include "util.c"
#ifndef ARGOPT_C__
#define ARGOPT_C__
#include "argopt.h"

/* initialize the argument structure */
argopt_t *argopt_open(void)
{
  argopt_t *ao;

  xnew(ao, 1);
  ao->narg = ao->nopt = 0;
  xnew(ao->args, 1);
  xnew(ao->opts, 1);
  return ao;
}

void argopt_close(argopt_t *ao)
{
  if (ao->args) free(ao->args);
  if (ao->opts) free(ao->opts);
  free(ao);
}

/* print help message and die */
static void argopt_help(argopt_t *ao)
{
  int i, len, maxlen;
  opt_t *lo = ao->opts;
  
  printf("%s [OPTIONS]", ao->prog);
  for (i = 0; i < ao->narg; i++) {
    if (strchr(ao->args[i].desc, ' ')) 
      printf(" (%s)", ao->args[i].desc);
    else
      printf(" %s", ao->args[i].desc);
  }
  printf("\n");
 
  printf("OPTIONS:\n") ;
  for (maxlen = 0, i = 0; i < ao->nopt; i++) {
    len =  strlen(lo[i].sflag);
    if (len > maxlen) maxlen = len;
  }
  for (i = 0; i < ao->nopt; i++)
    printf("%-*s : %s%s\n", maxlen, lo[i].sflag,
        (lo[i].hasarg ? "followed by " : ""), lo[i].desc);
  argopt_close(ao);
  exit(1);
}

/* register argument */
void argopt_regarg(argopt_t *ao, const char *fmt, void *ptr,
    const char *desc)
{
  arg_t *al = ao->args + ao->narg;

  al->fmt = fmt;
  al->ptr = ptr;
  al->desc = desc;
  if (fmt == NULL) *(const char **) ptr = NULL;
  ao->narg++;
  xrenew(ao->args, ao->narg+1);
}

/* register option */
void argopt_regopt(argopt_t *ao, const char *sflag, int hasarg,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *ol = ao->opts + ao->nopt;

  ol->sflag = sflag;
  ol->ch = sflag[2] ? '\0' : sflag[1]; /* no ch for a long flag */
  ol->hasarg = hasarg;
  ol->desc = desc;
  ol->fmt = fmt;
  ol->ptr = ptr;
  if (!hasarg) *(int *)ptr = 0;
  ao->nopt++;
  xrenew(ao->opts, ao->nopt+1);
}

static int sval_parse_(const char *val, const char *fmt, void *ptr)
{
  if (fmt == NULL) { /* raw string assignment */
    *(const char **)ptr = val;
  } else { /* call sscanf */
    int ret = sscanf(val, fmt, ptr);
    if (ret != 1) {
      fprintf(stderr, "unable to convert a value as fmt [%s], raw: [%s]\n",
          fmt, val);
      return 1;
    }
  }
  return 0;
}

/* main parser of arguments */
void argopt_parse(argopt_t *ao, int argc, const char **argv) 
{
  int i, j, k, ch, acnt = 0;
  arg_t *al = ao->args;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      if (acnt >= ao->narg) argopt_help(ao);
      al[acnt].val = argv[i];
      sval_parse_(al[acnt].val, al[acnt].fmt, al[acnt].ptr);
      ++acnt;
      continue;
    }

    /* it's an option, loop for abbreviated form "-abc" == "-a -b -c" */
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      int islong = (j == 1 && argv[i][1] == '-');
      
      if (islong) { /* compare against long options */
        for (k = 0; k < ao->nopt; k++)
          if (strncmp(argv[i], ol[k].sflag, strlen(ol[k].sflag)) == 0)
            break;
      } else { /* compare against short options */
        for (k = 0; k < ao->nopt; k++)
          if (ch == ol[k].ch)
            break;
      }
      if (k >= ao->nopt) {
        fprintf(stderr, "cannot handle option [%s]\n", argv[i]);
        argopt_help(ao);
      }
      if (ol[k].hasarg) { /* look for the additional argument for this */
        int hasv = 0;
        if (islong) { /* e.g., --version=11 */
          j = strlen(ol[k].sflag);
          if (argv[i][ j ] == '=') {
            ol[k].val = argv[i] + j + 1;
            hasv = 1;
          }
        } else { /* e.g., -n8 */
          if (argv[i][++j]) {
            ol[k].val = argv[i] + j;
            hasv = 1;
          }
        }
        
        if (!hasv) { /* --version 11 or -n 8 */
          if (++i >= argc) {
            printf("[%s|%s] requires an argument!\n", ol[k].sflag, argv[i-1]);
            argopt_help(ao);
          }
          ol[k].val = argv[i];
        }
        sval_parse_(ol[k].val, ol[k].fmt, ol[k].ptr);
        break; /* go to the next argument argv[i+1] */
      } else {
        ol[k].set = 1;
        *(int *)(ol[k].ptr) = 1;
        if (islong) break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
}

#endif
