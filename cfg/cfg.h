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

#endif

