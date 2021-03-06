#include "util.h"
#include "opt.h"
#include "ss.h"
#ifndef CFG_H__
#define CFG_H__



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



#define cfg_set(cfg, var) opt_isset(cfg->opts, cfg->nopt, &var, #var)


/* Read the value of a given variable from the current configuration file,
 * the name of variable is given by `key',
 * If the key is matched, its value is saved to `*var' through sscanf,
 *   otherwise, the content in *var is not modified.
 * If the function succeeds, it returns 0.
 * In case fmt is "%s", (*var) is a string, or a pointer to char.
 *   The space for (*var) will be managed through sscpy. */
INLINE int cfg_get(cfg_t *cfg, void *var, const char *key, const char *fmt)
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



/* load the whole configuration file into memory, parse it to entries */
INLINE cfg_t *cfg_open(const char *fn)
{
  cfg_t *cfg;
  cfgln_t *cln;
  FILE *fp;
  size_t i, j, n, size = 0;
  char *p, *q;

  xfopen(fp, fn, "r", return NULL);
  xnew(cfg, 1);
  if (ssfgetall(cfg->buf, &size, fp) == NULL) {
    fprintf(stderr, "error reading file %s\n", fn);
    return NULL;
  }
  sscat(cfg->buf, "\n"); /* in case the file is not ended by a new line, we add one */
  fclose(fp);

  /* count the number of lines (before allocating the key-table) */
  for (p = cfg->buf, i = 0, n = 0; i < size; i++) {
    if (p[i] == '\n' || p[i] == '\r') {
      if (i > 0 && p[i - 1] == '\\') {
        /* splice multiple lines: replace CR/LF with spaces
           the parser should be insensitive to the extra spaces */
        p[i - 1] = p[i] = ' ';
      } else { /* mark the end the current line */
        p[i] = '\0';
        n++;
      }

      /* replace successive CR/LF by spaces for efficiency
         the size of the key-table == the number of non-blank lines */
      for (j = i + 1; j < size && cisspace(p[j]); j++) p[j] = ' ';
    }
  }

  xnew(cfg->lns, n);

  /* load lines into the key-table */
  for (p = cfg->buf, j = 0, i = 0; i < size; i++) {
    if (cfg->buf[i] == '\0') {
      cfg->lns[j++].key = p;
      if (j >= n) break;
      p = cfg->buf + i + 1;
    }
  }
  n = j;

  /* parse each line to a key/value pair */
  for (j = 0; j < n; j++) {
    cln = cfg->lns + j;
    p = cln->key;
    strip(p); /* remove leading and trailing spaces */

    /* skip a blank or comment line */
    if (p[0] == '\0' || strchr("#;", p[0]) != NULL) {
      cln->key = NULL;
      continue;
    }

    /* remove trailing spaces and ';' */
    for ( q = p + strlen(p) - 1;
          q >= p && (cisspace(*q) || *q == ';'); q--)
      *q = '\0';

    /* skip a line without '=' */
    if ((q = strchr(p, '=')) == NULL) {
      cln->key = NULL;
      continue;
    }
    *q = '\0';
    cln->val = q + 1;
    strip(cln->key);
    strip(cln->val);
  }
  cfg->nln = (int) n;
  cfg->nopt = 0;
  xnew(cfg->opts, 1); /* s.t. we call `xrenew' the next time */
  return cfg;
}



INLINE void cfg_close(cfg_t *cfg)
{
  ssdelete(cfg->buf);
  free(cfg->lns);
  free(cfg->opts);
  memset(cfg, 0, sizeof(*cfg));
  free(cfg);
}



/* register an option request, return the index */
INLINE int cfg_add(cfg_t *cfg, const char *key, const char *fmt,
            void *ptr, const char *desc)
{
  /* if fmt == NULL, the memory space of (*ptr) will be invalid
   * after cfg_close(), "%s" is much safer */
  if (fmt == NULL) fmt = "%s";
  xrenew(cfg->opts, cfg->nopt + 1);
  opt_set(cfg->opts + cfg->nopt, NULL, key, fmt, ptr, desc);
  return cfg->nopt++;
}



#define CFG_UNUSED   0x0001
#define CFG_NOTSET   0x0002
#define CFG_CHECKUSE 0x0100
#define CFG_VERBOSE  0x1000
#define cfg_parse(cfg, flags) cfg_match(cfg, flags)

/* match requested options with entries in cfg file
 * returns 0 if successful
 * if mandatory variables are not set, the return `ret' contains CFG_NOTSET
 * if `flags' has OPT_CHECKUSE, the return `ret' has CFG_UNUSED if
 * there are unused variables */
INLINE int cfg_match(cfg_t *cfg, unsigned flags)
{
  int i, j, ret = 0, verbose = flags & CFG_VERBOSE, must;
  opt_t *o;
  cfgln_t *cln;

  for (i = 0; i < cfg->nopt; i++) {
    o = cfg->opts + i;
    for (j = 0; j < cfg->nln; j++) {
      cln = cfg->lns + j;
      if (cln->key != NULL && strcmp(cln->key, o->key) == 0) {
        cln->used = 1;
        o->flags |= OPT_SET;
        o->val = cln->val;
        opt_getval(o);
        break;
      }
    }
    must = (o->flags & OPT_MUST);
    if (!(o->flags & OPT_SET) && (must || verbose)) {
      fprintf(stderr, "cfg: %s not set, default: ", o->key);
      opt_fprintptr(stderr, o);
      fprintf(stderr, "\n");
      if (must) ret |= CFG_NOTSET;
    }
  }

  if (flags & CFG_CHECKUSE) {
    for (j = 0; j < cfg->nln; j++) {
      cln = cfg->lns + j;
      if (cln->key != NULL && !cln->used && verbose) {
        fprintf(stderr, "cfg: unused line: %s = %s\n",
            cln->key, cln->val);
        ret |= CFG_UNUSED;
      }
    }
  }
  return ret;
}



/* dump the current values */
INLINE void cfg_dump(const cfg_t *cfg)
{
  int i, len = 2;
  opt_t *ol = cfg->opts;

  /* get the width of the widest entry */
  for (i = 0; i < cfg->nopt; i++)
    len = intmax(len, strlen(ol[i].key));

  /* print values of all options */
  for (i = 0; i < cfg->nopt; i++) {
    fprintf(stderr, "%*s: ", len + 1, ol[i].key);
    opt_fprintptr(stderr, ol + i);
    fprintf(stderr, ",  %s\n", ol[i].desc);
  }
}

#endif

