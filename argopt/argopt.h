#include "util.h"
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



#define argopt_getopt(ao, p) opt_find(ao->opts, ao->nopt, p)
#define argopt_getarg argopt_getopt

/* test if argument/option is explicitly set */
#define argopt_isset(ao, var) opt_isset(ao->opts, ao->nopt, &var, #var)



/* initialize the argument structure */
INLINE argopt_t *argopt_open(unsigned flags)
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



INLINE void argopt_close(argopt_t *ao)
{
  if (ao->opts) { free(ao->opts); ao->opts = NULL; }
  free(ao);
}



/* print version and die */
INLINE void argopt_version(argopt_t *ao)
{
  fprintf(stderr, "%s: %s, version %d\n",
      ao->prog, ao->desc ? ao->desc : "", ao->version);
  if (ao->author && ao->tm)
    fprintf(stderr, "Copyright (c) %s %d\n", ao->author, ao->tm->tm_year + 1900);
  argopt_close(ao);
  exit(1);
}



/* print help message and die */
INLINE void argopt_help(argopt_t *ao)
{
  int i, len, maxlen;
  opt_t *o;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  fprintf(stderr, "%s, version %d",
      ao->desc ? ao->desc : ao->prog, ao->version);
  if (ao->author && ao->tm)
    fprintf(stderr, ", Copyright (c) %s %d", ao->author, ao->tm->tm_year + 1900);
  fprintf(stderr, "\nUSAGE\n  %s {OPTIONS}", ao->prog);
  for (i = 0; i < ao->nopt; i++) {
    const char *bra = "", *ket = "";
    o = ao->opts + i;
    if (o->isopt) continue;
    if (o->flags & OPT_MUST) {
      if (strchr(o->desc, ' '))
        bra = "[", ket = "]";
    } else
      bra = "{", ket = "}";
    fprintf(stderr, " %s%s%s", bra, o->desc, ket);
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "OPTIONS:\n") ;
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
    fprintf(stderr, "  %-*s : %s%s%s", maxlen, o->sflag,
        ((o->flags & OPT_MUST) ? "[MUST] " : ""),
        (!(o->flags & OPT_SWITCH) ? "followed by " : ""), desc);
    if (o->ptr && o->ptr != ao->dum_) { /* print default values */
      fprintf(stderr, ", default: ");
      opt_fprintptr(stderr, o);
    }
    fprintf(stderr, "\n");
  }
  argopt_close(ao);
  exit(1);
}



#define argopt_regarg(ao, fmt, ptr, desc) argopt_add(ao, NULL, fmt, ptr, desc)
#define argopt_regopt argopt_add
#define argopt_reghelp argopt_addhelp
#define argopt_regversion argopt_addversion
#define argopt_addhelp(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$HELP")
#define argopt_addversion(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$VERSION")

/* register an argument or option
 * sflag: string flag, or NULL for an argument
 * fmt: sscanf() format string, "%b" for a switch, "%r" for real
 * return the index */
INLINE int argopt_add(argopt_t *ao, const char *sflag,
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
INLINE void argopt_parse(argopt_t *ao, int argc, char **argv)
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

      if (islong) { /* match against long options */
        for (k = 0; k < ao->nopt; k++) {
          int lenf = strlen(ol[k].sflag);
          if (ol[k].isopt &&
              strncmp(argv[i], ol[k].sflag, lenf) == 0 &&
              strchr("= ", argv[i][lenf])) /* followed by a space or "=" */
            break;
        }
      } else { /* match against short options */
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

      if (ol[k].flags & OPT_SWITCH) { /* handle switches "%b" */
        ol[k].flags |= OPT_SET;
        /* switch the default value, note this flag may be passed
         * several times, so we don't want to flip around */
        *((int *) ol[k].ptr) = !ol[k].ival;
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
            fprintf(stderr, "%s(%s) requires an argument!\n", ol[k].sflag, argv[i - 1]);
            argopt_help(ao);
          }
          ol[k].val = argv[i];
        }
        ol[k].flags |= OPT_SET;
        if (0 != opt_getval(ol + k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & OPT_MUST) && !(ol[i].flags & OPT_SET)) {
      fprintf(stderr, "Error: missing %s %s: %s\n\n",
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

  /* get the width of the widest option */
  for (i = 0; i < ao->nopt; i++)
    if (ol[i].sflag)
      len = intmax(len, strlen(ol[i].sflag));

  /* print values of all options */
  for (i = 0; i < ao->nopt; i++) {
    const char *sflag = ol[i].sflag;
    if (sflag == NULL) sflag = "arg";
    fprintf(stderr, "%*s: ", len + 1, sflag);
    opt_fprintptr(stderr, ol + i);
    fprintf(stderr, ",  %s\n", ol[i].desc);
  }
}


#endif
