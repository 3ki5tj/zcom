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


unsigned cfgcheck(cfg_t *cfg, unsigned flags)
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


