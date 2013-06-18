#include "hist.h"

#define LCNT 3
#define XMIN 0.0
#define XMAX 1.0
#define XDEL 0.001

int main(void)
{
  int i, j, n = 10000000;
  unsigned wflags = HIST_VERBOSE | HIST_ADDAHALF; /* | HIST_KEEPHIST; */
  double x[LCNT] = {0.51, 0.51, 0.51}, lam[LCNT] = {3.8, 3.9, 3.999999};
  hist_t *hs;
  const char *fnhs = "hs.dat", *fnhs2 = "hs2.dat";

  hs = hs_open(LCNT, XMIN, XMAX, XDEL);
  /* generate histogram */
  for (i = 0; i < n; i++) {
    for (j = 0; j < LCNT; j++)
      x[j] = lam[j]*x[j]*(1 - x[j]);
    hs_add(hs, x, 1., HIST_VERBOSE);
  }
  hs_save(hs, fnhs, wflags);

  /* now try to load histogram */
  die_if (hs_load(hs, fnhs, HIST_VERBOSE | HIST_ADDITION) != 0,
    "cannot load histogram from %s\n", fnhs);
  /* write again */
  hs_save(hs, fnhs2, wflags);
  hs_close(hs);
  return 0;
}

