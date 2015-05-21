/* test on harmonic oscillator */
#include "md.h"
#include "av.h"
#include "hist.h"

int main(void)
{
  real x = 1, v = 1, f, dt = 0.002, zeta = 0, zeta2 = 1.f, tp = 1.;
  int t;
  hist_t *hs;

  hs = hist_open(1, -10.0, 10.0, 0.01);
  for (t = 1; t <= 1000000; t++) {
    v += f * dt * .5f;
    x += v * dt;
    f = -x;
    v += f * dt * .5f;
    md_vslang(&v, 1, 1, tp, dt, &zeta, zeta2, 1.f, NULL, NULL);
/*
    md_vrescalex(&v, 1, 1, tp, dt, NULL, NULL);
*/
    if (t % 100 == 0) {
      printf("%g %g\n", x, v);
    }
    hist_add1ez(hs, x, 0);
  }
  hist_save(hs, "x.his", HIST_ADDAHALF);
  hist_close(hs);
  return 0;
}
