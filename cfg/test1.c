#include "cfg.h"

int main(void)
{
  int   nr, cnt;
  float tmin, tmax;
  char  *p = NULL;
  cfg_t *cfg;
  const char *fncfg = "foo.cfg";

  die_if ((cfg = cfg_open(fncfg)) == NULL,
    "cannot read %s\n", fncfg);
  cfg_get(cfg, &nr,   "nrtemp", "%d");
  cfg_get(cfg, &tmin, "tmin",   "%f");
  cfg_get(cfg, &tmax, "tmax",   "%f");
  cfg_get(cfg, &p,    "mystr",  "%s");
  cfg_get(cfg, &cnt,  "arrcnt", "%d");

  printf("nr=%d, (%g,%g), cnt=%d\n", nr, tmin, tmax, cnt);
  printf("mystr=\"%s\"\n", p);
  cfg_close(cfg);
  ssdelall();
  return 0;
}

