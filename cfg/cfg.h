#include "util.h"
#include "opt.h"
#ifndef CFG_H__
#define CFG_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* a line in the cfg file: `key = val'
 * `key' and `val' point to spaces allocated in cfg->buf
 * and no additional memory are allocated */
typedef struct {
  char *key, *val;
  int used;
} cfgln_t;

/* the entire cfg file */
typedef struct {
  char *buf;      /* string buffer for the entire file */
  int nln;        /* number of input lines */
  cfgln_t *lns;   /* parsed lines */
  int nopt;       /* number of user-requested options */
  opt_t *opts;    /* user-requested options, such that
                     cfg works like command-line options */
} cfg_t;

typedef cfg_t cfgdata_t; /* old alias */

#define CFG_CHECKUSE 0x0100
#define CFG_VERBOSE  0x1000


cfg_t *cfg_open(const char *fn);
void cfg_close(cfg_t *cfg);
int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc);

#define CFG_UNUSED   0x0001
#define CFG_NOTSET   0x0002
#define cfg_parse(cfg, flags) cfg_match(cfg, flags)
int cfg_match(cfg_t *cfg, unsigned flags);

#define cfg_set(cfg, var) opt_isset(cfg->opts, cfg->nopt, &var, #var)

/* old style functions */
#define cfgopen(fn) cfg_open(fn)
#define cfgclose(cfg) cfg_close(cfg)
/* Read the value of a given variable from the current configuration file,
 * the name of variable is given by `key',
 * If the key is matched, its value is saved to `*var' through sscanf,
 *   otherwise, the content in *var is not modified.
 * If the function succeeds, it returns 0.
 * In case fmt is "%s", (*var) is a string, or a pointer to char.
 *   The space for (*var) will be managed through sscpy. */
INLINE int cfgget(cfg_t *cfg, void *var, const char *key, const char *fmt)
{
  int i;

  for (i = 0; i < cfg->nln; i++) {
    cfgln_t *cln = cfg->lns + i;
    if (cln->key != NULL && strcmp(cln->key, key) == 0) {
      if (strcmp(fmt, "%s") == 0) { /* string */
        sscpy( *((char **) var), cln->val); /* make a copy and return */
        return 0;
      } else /* use sscanf for other cases, like int, float,... */
        return EOF == sscanf(cln->val, fmt, var) ? 2 : 0;
    }
  }
  return 1; /* no match */
}
#endif  /* CFG_H__ */

