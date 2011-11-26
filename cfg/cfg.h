#include "util.h"
#include "opt.h"
#ifndef CFG_H__
#define CFG_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
  char *key, *val;
  int used;
} cfgent_t; /* line from cfg file */

typedef struct {
  char *buf;      /* the entire configuration file */
  int nent;       /* number of entries */
  cfgent_t *ents; /* entries */
  int nopt;       /* number of user-requested options */
  opt_t *opts;    /* user-requested options */
} cfg_t;
typedef cfg_t cfgdata_t;

#define CFG_CHECKDUP 0x0100
#define CFG_CHECKUSE 0x0200
#define CFG_VERBOSE  0x1000

#define cfgopen(fn) cfg_open(fn)
#define cfgclose(cfg) cfg_close(cfg)
#define cfgget(cfg, var, key, fmt) cfg_get(cfg, var, key, fmt)

cfg_t *cfg_open(const char *fn);
void cfg_close(cfg_t *cfg);
int cfg_get(cfg_t *cfg, void *var, const char *key, const char *fmt);
int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc);
int cfg_match(cfg_t *cfg, unsigned flags);

#endif  /* CFG_H__ */

