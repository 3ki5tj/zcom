#include "cfg.c"

int main(void)
{
  int   nr = 1, cnt = 0;
  float tmin = 0.5f, tmax = 1.0f;
  char  *p = NULL;
  cfg_t *cfg;

  if ((cfg = cfg_open("foo.cfg")) == NULL) {
    printf("error reading\n");
    return 1;
  }
  cfg_add(cfg, "nrtemp", "%d", &nr, "# of T.");
  cfg_add(cfg, "tmin", "%f", &tmin, "T min");
  cfg_add(cfg, "tmax", "%f", &tmax, "T max");
  cfg_add(cfg, "mystr",  "%s", &p, "some string");
  cfg_add(cfg, "arrcnt", "%d", &cnt, "array count");
  cfg_match(cfg, CFG_CHECKUSE);
  cfg_dump(cfg);

  printf("nr=%d, (%g,%g), cnt=%d\n", nr, tmin, tmax, cnt);
  printf("mystr=\"%s\"\n", p);
  cfg_close(cfg);
  ssdelall();
  return 0;
}

