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
  ao->nopt = 0;
  ao->opts = NULL;
  tmcmpl = time(NULL);
  ao->tm = localtime( &tmcmpl );
  memset(ao->dum_, '\0', sizeof(ao->dum_));
  return ao;
}

void argopt_close(argopt_t *ao)
{
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
  opt_t *o;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  printf("%s, version %d", 
      ao->desc ? ao->desc : ao->prog, ao->version);
  if (ao->author && ao->tm)
    printf(", Copyright (c) %s %d", ao->author, ao->tm->tm_year + 1900);
  printf("\nUSAGE\n  %s [OPTIONS]", ao->prog);
  for (i = 0; i < ao->nopt; i++) {
    const char *bra = "", *ket = "";
    o = ao->opts + i;
    if (o->isopt) continue;
    if (o->flags & OPT_MUST) {
      if (strchr(o->desc, ' ')) 
        bra = "{", ket = "}";
    } else
      bra = "[", ket = "]";
    printf(" %s%s%s", bra, o->desc, ket);
  }
  printf("\n");
 
  printf("OPTIONS:\n") ;
  for (maxlen = 0, i = 0; i < ao->nopt; i++) { /* compute the longest option */
    if (!ao->opts[i].isopt) continue;
    len = strlen(ao->opts[i].sflag);
    if (len > maxlen) maxlen = len;
  }
  for (i = 0; i < ao->nopt; i++) {
    o = ao->opts + i;
    if (!o->isopt) continue;
    desc = o->desc;
    if (strcmp(desc, "$HELP") == 0)
      desc = sysopt[0];
    else if (strcmp(desc, "$VERSION") == 0)
      desc = sysopt[1];
    printf("  %-*s : %s%s", maxlen, o->sflag,
        (!(o->flags & OPT_SWITCH) ? "followed by " : ""), desc);
    if (o->ptr && o->ptr != ao->dum_) { /* print default values */
      printf(", default: ");
      opt_printptr(o);
    }
    printf("\n");
  }
  argopt_close(ao);
  exit(1);
}

/* register an argument or option
 * sflag: string flag, or NULL for an argument
 * fmt: sscanf() format string, "%b" for a switch, "%r" for real
 * return the index */
int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *o;
  int n;

  n = ao->nopt++;
  xrenew(ao->opts, ao->nopt);
  o = ao->opts + n;
  opt_set(o, sflag, NULL, fmt, ptr, desc);
  return n;
}

/* main parser of arguments */
void argopt_parse(argopt_t *ao, int argc, char **argv) 
{
  int i, j, k, ch, acnt = 0;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      while (ol[acnt].isopt && acnt < ao->nopt) acnt++;
      if (acnt >= ao->nopt) argopt_help(ao);
      ol[acnt].val = argv[i];
      ol[acnt].flags |= OPT_SET;
      if (0 != opt_getval(ol + acnt))
        argopt_help(ao);
      ++acnt;
      continue;
    }

    /* it's an option, loop for abbreviated form "-abc" == "-a -b -c" */
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      int islong = (j == 1 && argv[i][1] == '-') | (ao->flags & ARGOPT_LONGOPT);
      
      if (islong) { /* compare against long options */
        for (k = 0; k < ao->nopt; k++)
          if (ol[k].isopt && 
              strncmp(argv[i], ol[k].sflag, strlen(ol[k].sflag)) == 0)
            break;
      } else { /* compare against short options */
        for (k = 0; k < ao->nopt; k++)
          if (ol[k].isopt && ch == ol[k].ch)
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
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & OPT_MUST) && !(ol[i].flags & OPT_SET)) {
      printf("Error: missing %s %s: %s\n\n",
          ol[i].isopt ? "option" : "argument", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}

/* dump the current values */
INLINE void argopt_dump(const argopt_t *ao)
{
  int i, len = 2;
  opt_t *ol = ao->opts;

  /* get the widest the option */
  for (i = 0; i < ao->nopt; i++)
    len = intmax(len, strlen(ol[i].sflag));

  /* print values of all options */
  for (i = 0; i < ao->nopt; i++) {
    printf("%*s: ", len+1, ol[i].sflag);
    opt_printptr(ol + i);
    printf(",  %s\n", ol[i].desc);
  }
}


#endif
