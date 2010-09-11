#ifndef CFG_H__
#define CFG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
  int    n;           /* number of lines */
  int    canfree;     /* whether the struct is dynamically allocated */
  char **key,**value; /* key[i] = value[i] */
  char  *buf;         /* the whole configuration file */
} cfgdata_t;

cfgdata_t *cfgopen(const char *filenm);
void cfgclose(cfgdata_t *cfg);
int cfgget(cfgdata_t *cfg, void *var, const char *key, const char* fmt);

#endif  /* CFG_H__ */




/* print a mismatch message and die */
#define CFG_MISMATCH_TPFMT_(fmt, tp) \
    fprintf(stderr, "format %s doesn\'t match type %s", fmt, #tp);    \
    exit(1);

/* check if a type matches format by sizeof */
#define CFG_CHECKTPFMT(tp, fmt)                                       \
  if (fmt[0] != '%') {                                                \
    fprintf(stderr, "invalid format string %s\n", fmt);               \
    exit(1);                                                          \
  } else if (strchr("di", fmt[1]) && strcmp(#tp, "int") != 0) {       \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("uxo", fmt[1]) && strcmp(#tp, "unsigned") != 0) { \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("efg", fmt[1]) && strcmp(#tp, "float") != 0) {    \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strcmp(fmt, "%lf") == 0 && strcmp(#tp, "double") != 0) { \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("sc", fmt[1]) && strcmp(#tp, "char *") != 0) {    \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  }

/* read variable var in format fmt */
#define CFG_GETATOM0_(var, name, tp, fmt) {       \
  CFG_CHECKTPFMT(tp, fmt)                         \
  err = cfgget(cfg, &(var), name, fmt); }


