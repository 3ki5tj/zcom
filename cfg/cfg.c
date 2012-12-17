#include "ss.c"
#ifndef CFG_C__
#define CFG_C__
#include "cfg.h"

/* load the whole configuration file into memory, parse it to entries */
cfg_t *cfg_open(const char *fn)
{
  cfg_t *cfg;
  cfgent_t *ent;
  FILE *fp;
  size_t i, j, n, size = 0;
  char *p, *q;

  xnew(cfg, 1);

  xfopen(fp, fn, "r", return NULL);
  if (ssfgetall(cfg->buf, &size, fp) == NULL) {
    fprintf(stderr, "error reading file %s\n", fn);
    return NULL;
  }
  sscat(cfg->buf, "\n"); /* in case the file is not ended by a new line, we add one */
  fclose(fp);

  /* count the number of lines (before allocating the key-table) */
  for (p = cfg->buf, i = 0, n = 0; i < size; i++) {
    if (p[i] == '\n' || p[i] == '\r') {
      if (i > 0 && p[i-1] == '\\') {
        /* multiple-line splicing by replacing cr, lf with spaces
           parse should be aware of these additional spaces */
        p[i-1] = p[i] = ' ';
      } else {
        p[i] = '\0';
        n++;
      }

      /* replace following CR LF by spaces for efficiency
         as the size of the key table == the number of blank lines */
      for (j = i+1; j < size && cisspace(p[j]); j++) p[j] = ' ';
    }
  }

  xnew(cfg->ents, n);

  /* load lines into the keytable */
  for (p = cfg->buf, j = 0, i = 0; i < size; i++) {
    if (cfg->buf[i] == '\0') {
      cfg->ents[j++].key = p;
      if (j >= n) break;
      p = cfg->buf + i + 1;
    }
  }
  n = j;

  /* parse each line to a key-value pair */
  for (j = 0; j < n; j++) {
    ent = cfg->ents + j;
    p = ent->key;
    strip(p);

    /* skip a blank or comment line */
    if (p[0] == '\0' || strchr("#%!;", p[0]) != NULL) {
      ent->key = NULL;
      continue;
    }

    /* remove trailing space and ';' */
    for (q = p + strlen(p) - 1; q >= p && (cisspace(*q) || *q == ';'); q--)
      *q = '\0';

    if ((q = strchr(p, '=')) == NULL) { /* skip a line without '=' */
      ent->key = NULL;
      continue;
    }
    *q = '\0';
    ent->val = q + 1;
    strip(ent->key);
    strip(ent->val);
  }
  cfg->nent = (int) n;
  cfg->nopt = 0;
  xnew(cfg->opts, 1);
  return cfg;
}

void cfg_close(cfg_t *cfg)
{
  ssdelete(cfg->buf);
  free(cfg->ents);
  free(cfg->opts);
  memset(cfg, 0, sizeof(*cfg));
  free(cfg);
}

/* register an option request, return the index */
int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc)
{
  int n = cfg->nopt++;
  opt_t *o;
  xrenew(cfg->opts, cfg->nopt);
  o = cfg->opts + n;
  opt_set(o, NULL, key, fmt, ptr, desc);
  return n;
}

/* match requested options with entries in cfg file
 * returns 0 if successful
 * if mandetory variables are not set  
 * if `flags' has OPT_CHECKUSE, then unused setting cause an error code CFG_UNUSED */
int cfg_match(cfg_t *cfg, unsigned flags)
{
  int i, j, ret = 0, verbose = flags & CFG_VERBOSE, must;
  opt_t *o;
  cfgent_t *ent;

  for (i = 0; i < cfg->nopt; i++) {
    o = cfg->opts + i;
    for (j = 0; j < cfg->nent; j++) {
      ent = cfg->ents + j;
      if (ent->key != NULL && strcmp(ent->key, o->key) == 0) {
        ent->used = 1;
        o->flags |= OPT_SET;
        o->val = ent->val;
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
    for (j = 0; j < cfg->nent; j++) {
      ent = cfg->ents + j;
      if (ent->key != NULL && !ent->used && verbose) {
        fprintf(stderr, "cfg: unused entry: %s = %s\n", ent->key, ent->val);
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
    fprintf(stderr, "%*s: ", len+1, ol[i].key);
    opt_fprintptr(stderr, ol + i);
    fprintf(stderr, ",  %s\n", ol[i].desc);
  }
}

#endif

