#include "util.c"
#ifndef ARGOPT_C__
#define ARGOPT_C__
#include "argopt.h"

/* initialize the argument structure */
argopt_t *argopt_open(int version, const char *desc, const char *author)
{
  argopt_t *ao;
  time_t tmcmpl;

  xnew(ao, 1);
  ao->narg = ao->nopt = 0;
  ao->version = version;
  ao->desc = desc;
  ao->author = author;
  tmcmpl = time(NULL);
  ao->tm = localtime( &tmcmpl );
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

/* print version and die */
static void argopt_version(argopt_t *ao)
{
  printf("%s: %s, version %d. Copyright (c) %s %d\n", 
      ao->prog, ao->desc ? ao->desc : "", 
      ao->version, (ao->author ? ao->author : ""), 
      ao->tm->tm_year + 1900);
  argopt_close(ao);
  exit(1);
}

/* print help message and die */
static void argopt_help(argopt_t *ao)
{
  int i, len, maxlen;
  opt_t *lo = ao->opts;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  printf("%s, version %d. Copyright (c) %s %d\n", 
      ao->desc ? ao->desc : ao->prog, 
      ao->version, (ao->author ? ao->author : ""), 
      ao->tm->tm_year + 1900);
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
  for (i = 0; i < ao->nopt; i++) {
    desc = lo[i].desc;
    if (strcmp(desc, "$HELP") == 0)
      desc = sysopt[0];
    else if (strcmp(desc, "$VERSION") == 0)
      desc = sysopt[1];
    printf("  %-*s : %s%s\n", maxlen, lo[i].sflag,
        (lo[i].hasarg ? "followed by " : ""), desc);
  }
  argopt_close(ao);
  exit(1);
}

/* register argument */
opt_t *argopt_regarg(argopt_t *ao, int must,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *al = ao->args + ao->narg;

  al->isopt = 0;
  al->ch = '\0';
  al->sflag = NULL;
  al->hasarg = 1;
  al->must = must;
  al->hasarg = 0;
  al->set = 0;
  al->fmt = fmt;
  if (ptr == NULL) {
    printf("using null pointer for argument %s\n", desc);
    exit(1);
  }
  al->ptr = ptr;
  al->desc = desc;
  
  ao->narg++;
  xrenew(ao->args, ao->narg+1);
  return al;
}

/* register option */
opt_t *argopt_regopt(argopt_t *ao, const char *sflag, int hasarg,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *ol = ao->opts + ao->nopt;

  ol->isopt = 1;
  ol->sflag = sflag;
  ol->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  ol->hasarg = hasarg;
  ol->must = 0; /* TODO: mandatory options */
  ol->set = 0;
  ol->desc = desc;
  ol->fmt = fmt;
  ol->ptr = ptr;
  if (ptr == NULL) {
    printf("using null pointer for option %s: %s\n", sflag, desc);
    exit(1);
  }
  ao->nopt++;
  xrenew(ao->opts, ao->nopt+1);
  return ol;
}

/* translate string values to actual ones through sscanf() */
static int opt_getval_(opt_t *o)
{
  const char *fmt = o->fmt;
  if (fmt == NULL) { /* raw string assignment */
    *(const char **)o->ptr = o->val;
  } else { /* call sscanf */
    int ret, ch;
    ch = fmt[strlen(fmt)-1];
    if (ch == 'r') /* real */
      fmt = (sizeof(real) == 4) ? "%f" : "%lf";
    ret = sscanf(o->val, fmt, o->ptr);
    if (ret != 1) {
      fprintf(stderr, "Error: unable to convert a value for [%s] as fmt [%s], raw string: [%s]\n",
          o->desc, fmt, o->val);
      return 1;
    }
  }
  return 0;
}

/* main parser of arguments */
void argopt_parse(argopt_t *ao, int argc, const char **argv) 
{
  int i, j, k, ch, acnt = 0;
  opt_t *al = ao->args;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      if (acnt >= ao->narg) argopt_help(ao);
      al[acnt].val = argv[i];
      al[acnt].set = 1;
      if (0 != opt_getval_(al+acnt))
        argopt_help(ao);
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
      
      if (ol[k].desc[0] == '$') { /* system commands */
        if (strcmp(ol[k].desc, "$HELP") == 0)
          argopt_help(ao);
        else if (strcmp(ol[k].desc, "$VERSION") == 0)
          argopt_version(ao);
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
        if (0 != opt_getval_(ol+k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      } else {
        ol[k].set = 1;
        *(int *)(ol[k].ptr) = 1;
        if (islong) break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->narg; i++) {
    if (al[i].must && !al[i].set) {
      printf("Error: missing argument %d: %s\n\n", i, al[i].desc);
      argopt_help(ao);
    }
  }
  for (i = 0; i < ao->nopt; i++) {
    if (ol[i].must && !ol[i].set) {
      printf("Error: missing option %s: %s\n\n", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}

#endif
