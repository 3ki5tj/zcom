#include "util.h"
#ifndef ARGOPT_C__
#define ARGOPT_C__
#include "argopt.h"

/* initialize the argument structure */
argopt_t *argopt_open(unsigned flags)
{
  argopt_t *ao;
  time_t tmcmpl;

  xnew(ao, 1);
  ao->flags = flags;
  ao->narg = ao->nopt = 0;
  ao->args = ao->opts = NULL;
  tmcmpl = time(NULL);
  ao->tm = localtime( &tmcmpl );
  memset(ao->dum_, '\0', sizeof(ao->dum_));
  return ao;
}

void argopt_close(argopt_t *ao)
{
  if (ao->args) { free(ao->args); ao->args = NULL; }
  if (ao->opts) { free(ao->opts); ao->opts = NULL; }
  free(ao);
}

/* print version and die */
static void argopt_version(argopt_t *ao)
{
  printf("%s: %s, version %d\n", 
      ao->prog, ao->desc ? ao->desc : "", ao->version);
  if (ao->author && ao->tm)
    printf("Copyright (c) %s %d\n", ao->author, ao->tm->tm_year + 1900);
  argopt_close(ao);
  exit(1);
}

/* print help message and die */
static void argopt_help(argopt_t *ao)
{
  int i, len, maxlen;
  opt_t *lo;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  printf("%s, version %d", 
      ao->desc ? ao->desc : ao->prog, ao->version);
  if (ao->author && ao->tm)
    printf(", Copyright (c) %s %d", ao->author, ao->tm->tm_year + 1900);
  printf("\nUSAGE\n  %s [OPTIONS]", ao->prog);
  for (i = 0; i < ao->narg; i++) {
    const char *bra = "", *ket = "";
    lo = ao->args + i;
    if (lo->flags & OPT_MUST) {
      if (strchr(lo->desc, ' ')) 
        bra = "{", ket = "}";
    } else
      bra = "[", ket = "]";
    printf(" %s%s%s", bra, lo->desc, ket);
  }
  printf("\n");
 
  printf("OPTIONS:\n") ;
  for (maxlen = 0, i = 0; i < ao->nopt; i++) {
    len =  strlen(ao->opts[i].sflag);
    if (len > maxlen) maxlen = len;
  }
  for (i = 0; i < ao->nopt; i++) {
    lo = ao->opts + i;
    desc = lo->desc;
    if (strcmp(desc, "$HELP") == 0)
      desc = sysopt[0];
    else if (strcmp(desc, "$VERSION") == 0)
      desc = sysopt[1];
    printf("  %-*s : %s%s", maxlen, lo->sflag,
        (!(lo->flags & OPT_SWITCH) ? "followed by " : ""), desc);
    if (lo->ptr && lo->ptr != ao->dum_) { /* print default values */
      printf(", default: ");
      opt_printptr(lo);
    }
    printf("\n");
  }
  argopt_close(ao);
  exit(1);
}

/* register option: fmt = "%b" for a switch
 * return the index */
int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *o;
  int n;

  if (sflag) { /* option */
    n = ao->nopt++;
    xrenew(ao->opts, ao->nopt);
    o = ao->opts + n;
  } else { /* argument */
    n = ao->narg++;
    xrenew(ao->args, ao->narg);
    o = ao->args + n;
  }
  opt_set(o, sflag, NULL, fmt, ptr, desc);
  return n;
}

/* main parser of arguments */
void argopt_parse(argopt_t *ao, int argc, char **argv) 
{
  int i, j, k, ch, acnt = 0;
  opt_t *al = ao->args;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      if (acnt >= ao->narg) argopt_help(ao);
      al[acnt].val = argv[i];
      al[acnt].flags |= OPT_SET;
      if (0 != opt_getval(al+acnt))
        argopt_help(ao);
      ++acnt;
      continue;
    }

    /* it's an option, loop for abbreviated form "-abc" == "-a -b -c" */
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      int islong = (j == 1 && argv[i][1] == '-') | (ao->flags & ARGOPT_LONGOPT);
      
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

      if (ol[k].flags & OPT_SWITCH) {
        ol[k].flags |= OPT_SET;
        *(int *)ol[k].ptr = 1;
        if (islong) break; /* go to the next argument argv[i+1] */
      } else { /* look for the additional argument for this */
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
            printf("%s(%s) requires an argument!\n", ol[k].sflag, argv[i-1]);
            argopt_help(ao);
          }
          ol[k].val = argv[i];
        }
        ol[k].flags |= OPT_SET;
        if (0 != opt_getval(ol+k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->narg; i++) {
    if ((al[i].flags & OPT_MUST) && !(al[i].flags & OPT_SET)) {
      printf("Error: missing argument %d: %s\n\n", i, al[i].desc);
      argopt_help(ao);
    }
  }
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & OPT_MUST) && !(ol[i].flags & OPT_SET)) {
      printf("Error: missing option %s: %s\n\n", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}

#endif
