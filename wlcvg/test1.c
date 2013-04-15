#include <stdlib.h>
#include "wlcvg.c"

int main(void)
{
  wlcvg_t *wl;
  int t, flags = WLCVG_UPDLNFC;
  double x = 0.5, x1;
  double percut = 0.2, lnfc = 0.0;

  wl = wlcvg_open(lnfc, 1.0, 0.5, percut, 0.0, 0, 1, 0.01, flags);
  for (t = 1; t <= 1000000; t++) {
    x1 = x + 0.02 * (2.*rand()/RAND_MAX - 1);
    if (x1 <= 1 && x1 >= 0) x = x1;
    wlcvg_update(wl, x);
    if (t % 10000 == 0)
      printf("t %d: stage %d, lnf %g, per %g\n", 
          t, wl->stage, wl->lnf, wl->perc);
  }
  wlcvg_close(wl);
  return 0;
}
