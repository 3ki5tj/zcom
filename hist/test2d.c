#include "hist2.h"

#define LCNT 3

#define XMIN 0.0
#define XMAX 1.0
#define XDEL 0.01

int main(void)
{
  int i, j, n = 1000000;
  unsigned wflags = HIST2_VERBOSE | HIST2_ADDAHALF | HIST2_KEEPHIST;
  double x[LCNT] = {0.51, 0.51, 0.51},
         y[LCNT] = {0., 0., 0.};
  hist2_t *hs;
  const char *fnhs = "hs.dat";

  hs = hs2_open(LCNT, XMIN, XMAX, XDEL, XMIN, XMAX, XDEL);
  /* generate histogram */
  for (i = 0; i < n; i++) {
    for (j = 0; j < LCNT; j++) {
      y[j] = 1. * rand() / RAND_MAX;
      x[j] = 1. * rand() / RAND_MAX;
    }
    hs2_add(hs, x, y, 1, 1., HIST2_VERBOSE);
  }
  hs2_save(hs, fnhs, wflags);

  /* now try to load histogram */
  if (0 != hs2_load(hs, fnhs, HIST2_VERBOSE)) {
    fprintf(stderr, "cannot load histogram\n");
    return -1;
  }
  /* write again */
  hs2_save(hs, fnhs, wflags);
  hs2_close(hs);
  return 0;
}

