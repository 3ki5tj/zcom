#include "cfg.c"

int main(void)
{
  int   nr, cnt;
  float tmin, tmax;
  char  *p = NULL;
  cfgdata_t *cfg;
  
  if ((cfg = cfgopen("foo.cfg")) == NULL) {
    printf("error reading\n");
    return 1;
  }
  cfgget(cfg, &nr,   "nrtemp", "%d");
  cfgget(cfg, &tmin, "tmin",   "%f");
  cfgget(cfg, &tmax, "tmax",   "%f");
  cfgget(cfg, &p,    "scode",  "%s");
  cfgget(cfg, &cnt,  "arrcnt", "%d");

  printf("nr=%d, (%g,%g), cnt=%d\n", nr, tmin, tmax, cnt);
  printf("scode=\"%s\"\n", p);
  cfgcheck(cfg, CFG_VERBOSE|CFG_CHECKDUP|CFG_CHECKUSE);
  cfgclose(cfg);
  ssdelall();
  return 0;
}

