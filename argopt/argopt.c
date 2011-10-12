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
  const char *sysopt[2] = {"print help message", "print version"}, *desc, *fmt;

  printf("%s, version %d", 
      ao->desc ? ao->desc : ao->prog, ao->version);
  if (ao->author && ao->tm)
    printf(", Copyright (c) %s %d", ao->author, ao->tm->tm_year + 1900);
  printf("\nUSAGE\n  %s [OPTIONS]", ao->prog);
  for (i = 0; i < ao->narg; i++) {
    const char *bra = "", *ket = "";
    lo = ao->args + i;
    if (lo->flags & ARGOPT_MUST) {
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
        (!(lo->flags & ARGOPT_SWITCH) ? "followed by " : ""), desc);
    if (lo->ptr && lo->ptr != &ao->dum_) { /* print default values */
      printf(", default: ");
      for (fmt = lo->fmt; *fmt && *fmt != '%'; fmt++) ;
#define ELIF_PF_(fm, fmp, type) else if (strcmp(fmt, fm) == 0) printf((lo->pfmt ? lo->pfmt : fmp), *(type *)lo->ptr)
      if (fmt == NULL || *fmt == '\0') printf("%s", (*(char **)lo->ptr) ? (*(char **)lo->ptr) : "NULL");
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
      else printf("unknown %s-->%%d: %d\n", fmt, *(int *)lo->ptr);
#undef ELIF_PF_
    }
    printf("\n");
  }
  argopt_close(ao);
  exit(1);
}

/* register option: fmt = "%b" for a switch */
int argopt_regopt(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *ol;
  int n;

  if (sflag) { /* option */
    n = ao->nopt++;
    xrenew(ao->opts, ao->nopt);
    ol = ao->opts + n;
    ol->isopt = 1;
    ol->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    n = ao->narg++;
    xrenew(ao->args, ao->narg);
    ol = ao->args + n;
    ol->isopt = 0;
    ol->ch = '\0';
  }
  ol->sflag = sflag;
  ol->flags = 0;
  die_if (ptr == NULL, "null pass to argopt with %s: %s\n", sflag, desc);
  ol->ptr = ptr;
  if (fmt == NULL) fmt = "";
  if (fmt[0] == '!') {
    fmt++;
    ol->flags |= ARGOPT_MUST;
  }
  if (strcmp(fmt, "%b") == 0) {
    fmt = "%d";
    ol->flags |= ARGOPT_SWITCH;
  }
  ol->fmt = fmt;
  ol->pfmt = NULL;
  ol->desc = desc;
  return n;
}

/* translate string values to actual ones through sscanf() */
static int opt_getval_(opt_t *o)
{
  const char *fmt = o->fmt;
  
  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
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
      al[acnt].flags |= ARGOPT_SET;
      if (0 != opt_getval_(al+acnt))
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

      if (ol[k].flags & ARGOPT_SWITCH) {
        ol[k].flags |= ARGOPT_SET;
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
        ol[k].flags |= ARGOPT_SET;
        if (0 != opt_getval_(ol+k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->narg; i++) {
    if ((al[i].flags & ARGOPT_MUST) && !(al[i].flags & ARGOPT_SET)) {
      printf("Error: missing argument %d: %s\n\n", i, al[i].desc);
      argopt_help(ao);
    }
  }
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & ARGOPT_MUST) && !(ol[i].flags & ARGOPT_SET)) {
      printf("Error: missing option %s: %s\n\n", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}

#endif
