#include "util.h"
#include "ss.c"

#ifndef CFG_C__
#define CFG_C__
/*
 * =====================================================================
 *
 * simple configuration file
 *
 * =====================================================================
 *
 * to open a configuration file:
 *    cfg = cfgopen("your.cfg"); 
 * to load a parameter,
 *    cfgget(cfg, &var, "var_name", scanf_fmt);
 * to finish up:
 *    cfgclose(cfg); 
 */
#include "cfg.h"

#define isspace_(c) isspace((unsigned char)(c))

/* load the whole configuration file into memory (return value);
 * parse it to entries, return 0 if successful */
static int cfgload_(cfgdata_t *cfg, const char *filenm)
{
  FILE *fp;
  size_t i, j, n, size = 0;
  char *p, *q, *lin;

  xfopen(fp, filenm, "r", return -1);
  if (ssfgetall(cfg->buf, &size, fp) == NULL) {
    fprintf(stderr, "error reading file %s\n", filenm);
    return 4;
  }
  sscat(cfg->buf, "\n"); /* in case the file is not ended by a new line, we add one */
  fclose(fp);

#ifdef CFGDBG_
  printf("the cfg is loaded, size=%d\n", size);
#endif

  /* count the number of lines (before allocating the key-table) */
  for (i = 0, n = 0; i < size; i++) {
    if (cfg->buf[i] == '\n' || cfg->buf[i] == '\r') {
      if (i > 0 && cfg->buf[i-1] == '\\') {
        /* multiple-line splicing by replacing cr, lf with spaces
           parse should be aware of these additional spaces */
        cfg->buf[i-1] = cfg->buf[i] = ' ';
      } else {
        cfg->buf[i] = '\0';
        n++;
      }

      /* we replace following cr & lf's by spaces for efficiency
         as the size of the key table == the number of blank lines */
      for (j = i+1; j < size && isspace_(cfg->buf[j]); j++)
        cfg->buf[j] = ' ';
#ifdef CFGDBG_
      if (j-1 >= i+1) printf("j=%d to %d are replaced by spaces\n", i+1, j-1);
#endif
    }
  }
#ifdef CFGDBG_
  printf("# of lines: %d\n", n);
#endif

  xnew(cfg->key, n);
  xnew(cfg->val, n);
  xnew(cfg->used, n);

  /* load lines into the keytable */
  for (p = cfg->buf, j = 0, i = 0; i < size; i++) {
    if (cfg->buf[i] == '\0') {
      cfg->key[j++] = p;
      if (j >= n) break;
      p = cfg->buf + i + 1;
    }
  }
  n = j;
#ifdef CFGDBG_
  fprintf(stderr, "load %d lines.\n", n);
#endif

  /* now parse lines: separate values from keys */
  for (j = 0; j < n; j++) {
    lin = cfg->key[j];

    /* remove the the leading spaces */
    for (; *lin && isspace_(*lin); lin++) ;
    cfg->key[j] = lin;
    /* skip a blank or comment line */
    if (lin[0] == '\0' || strchr("#%!;", lin[0]) != NULL) {
      cfg->key[j] = NULL;
      continue;
    }

    /* remove trailing space and ';' */
    for (q = lin + strlen(lin) - 1;
         q >= lin && (isspace_(*q) || *q == ';'); q--) *q = '\0';

    if ((q = strchr(lin, '=')) == NULL) { /* skip a line without '=' */
      cfg->key[j] = NULL;
      continue;
    }

    /* find the end of key --> 'q' */
    *q = '\0';
    p  = q + 1;
    for (--q; isspace_(*q); q--) *q = '\0';
    for (; (*p) && isspace_(*p); p++) ; /* skip leading space, 'p' -> value */
    cfg->val[j] = p;
  }

#ifdef CFGDBG_
  for (j = 0; j < n; j++) {
    if (cfg->key[j] != NULL)
      printf("key=%s, value=%s\n", cfg->key[j], cfg->val[j]);
  }
  printf("%d lines\n", n);
  getchar();
#endif

  cfg->n = (int) n;
  return 0;
}

#undef isspace_

/* a wrapper of cfgload_ to make it more like fopen */
cfgdata_t *cfgopen(const char *filenm)
{
  cfgdata_t *cfg;

  if ((cfg = calloc(1, sizeof(*cfg))) == NULL) {
    fprintf(stderr, "cannot allocate space for cfgdata_t.\n");
    return NULL;
  }
  if (cfgload_(cfg, filenm) != 0) {
    free(cfg);
    return NULL;
  }
  cfg->canfree = 1; /* so it can be safely freed */
  return cfg;
}


void cfgclose(cfgdata_t *cfg)
{
  int canfree = cfg->canfree; /* save the value before memset */

  free(cfg->used);
  free(cfg->val);
  free(cfg->key);
  ssdelete(cfg->buf);
  memset(cfg, 0, sizeof(*cfg));
  if (canfree) free(cfg); /* free cfg if it is created by cfgopen */
}


/* Read the value of a given variable from the current configuration file,
 * the name of variable is given by `key',
 * If the key is matched, its value is saved to `*var' through sscanf,
 *   otherwise, the content in *var is not modified.
 * If the function succeeds, it returns 0.
 * A comment line in the configuration file starts with '#', '%' or '!'.
 * In case fmt is "%s", (*var) is a string, or a pointer to char.
 *   The space for (*var) will be managed through sscpy, which should *not*
 *   to be free'd.  Instead, ssdel should be called if really necessary.
 * */
int cfgget(cfgdata_t *cfg, void *var, const char *key, const char* fmt)
{
  int j;

  if (cfg == NULL || cfg->key == NULL || var == NULL || key == NULL || fmt == NULL) {
    fprintf(stderr, "cfgget: NULL pointer.\n");
    return 1;
  }

  for (j = 0; j < cfg->n; j++)
    if (cfg->key[j] != NULL && strcmp(cfg->key[j], key) == 0) {
      cfg->used[j] = 1;
      if (strcmp(fmt, "%s") == 0) { /* string case */
        sscpy( *(char **)var, cfg->val[j]); /* make a copy and return */
        return 0;
      } else { /* use sscanf for other cases, like int, float,... */
        return (EOF == sscanf(cfg->val[j], fmt, var)) ? 2 : 0;
      }
    }
  return 1; /* no match */
}

unsigned cfgcheck(cfgdata_t *cfg, unsigned flags)
{
  int i, j, n = cfg->n;
  unsigned verbose = (flags & CFG_VERBOSE), ret = 0;

  if (flags & CFG_CHECKDUP) { /* check duplicated keys */
    for (i = 0; i < n; i++) {
      if (cfg->key[i] == NULL) continue;
      for (j = i+1; j < n; j++) {
        if (cfg->key[j] == NULL) continue;
        if (strcmp(cfg->key[i], cfg->key[j]) == 0) {
          if (verbose)
            fprintf(stderr, "duplicated setup, key %s: %s vs %s\n",
                cfg->key[i], cfg->val[i], cfg->val[j]);
          ret |= CFG_CHECKDUP;
        }
      }
    }
  }
  
  if (flags & CFG_CHECKUSE) {
    for (i = 0; i < n; i++) {
      if (cfg->key[i] == NULL || cfg->used[i]) continue;
      if (verbose)
        fprintf(stderr, "unused key %s: %s\n", cfg->key[i], cfg->val[i]);
      ret |= CFG_CHECKUSE;
    }
  }
  return ret;
}

#ifdef CFG_LEGACY
/* read the value of a given array variable from the current configuration file,
   the size of each array item is given by `itemsize',
   the name of the array starts with `key' followed by an index
   e.g. if key="arr", entries in configuration files are `arr0', `arr1', ...
   indices are from i0 to iN-1;
   if the function succeeds, it returns 0. */
int cfggetarr(cfgdata_t *cfg, void const *varr, size_t itemsize,
    const char *key, const char* fmt, int i0, int iN)
{
  int i;
  char *var, itemname[128];
  var = (char*)varr;
  if (strlen(key) > sizeof(itemname)-16) {
    fprintf(stderr, "key name is too long\n");
    return 1;
  }
  for (i = i0; i < iN; i++) {
    sprintf(itemname, "%s%d", key, i);
    cfgget(cfg, var, itemname, fmt);
    var += itemsize;
  }
  return 0;
}
#endif

#endif

