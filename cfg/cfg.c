#include "ss.c"
#ifndef CFG_C__
#define CFG_C__
#include "cfg.h"

/* load the whole configuration file into memory;
 * parse it to entries */
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
  free(cfg); /* free cfg if it is created by cfgopen */
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
int cfg_get(cfg_t *cfg, void *var, const char *key, const char *fmt)
{
  int i;

  for (i = 0; i < cfg->nent; i++) {
    cfgent_t *ent = cfg->ents + i;
    if (ent->key != NULL && strcmp(ent->key, key) == 0) {
      if (strcmp(fmt, "%s") == 0) { /* string case */
        sscpy( *(char **)var, ent->val); /* make a copy and return */
        return 0;
      } else { /* use sscanf for other cases, like int, float,... */
        return (EOF == sscanf(ent->val, fmt, var)) ? 2 : 0;
      }
    }
  }
  return 1; /* no match */
}

/* add option, return the index */
int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc)
{
  opt_t *o;
  int n = cfg->nopt++;
  xrenew(cfg->opts, cfg->nopt);
  o = cfg->opts + n;
  opt_set(o, NULL, key, fmt, ptr, desc);
  return n;
}

/* match requested options with cfg file entries */
int cfg_match(cfg_t *cfg, unsigned flags)
{
  int i, j, ret = 0;
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
    if (!(o->flags & OPT_SET)) {
      printf("cfg: %s not set, default: ", o->key);
      opt_printptr(o);
      printf("\n");
      if (o->flags & OPT_MUST) ret = 1;
    }
  }

  if (flags & CFG_CHECKUSE) {
    for (j = 0; j < cfg->nent; j++) {
      ent = cfg->ents + j;
      if (ent->key != NULL && !ent->used)
        printf("cfg: unused entry: %s = %s\n", ent->key, ent->val);
    }
  }
  return ret;
}

#endif

