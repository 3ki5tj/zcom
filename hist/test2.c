#include "hist.c"

#define LCNT 3

#define XMIN 0.0
#define XMAX 1.0
#define XDEL 0.01

int main(void)
{
  int i, j, n = 1000000;
  unsigned wflags = HIST_VERBOSE | HIST_ADDAHALF | HIST_KEEPHIST;
  double x[LCNT] = {0.51, 0.51, 0.51},
         y[LCNT] = {0., 0., 0.};
  //double lam[LCNT] = {3.8, 3.9, 3.999999};
  hist2_t *hs;

  hs = hs2_opensqr(LCNT, XMIN, XMAX, XDEL);
  /* generate histogram */
  for (i = 0; i < n; i++) {
    for (j = 0; j < LCNT; j++) {
      y[j] = 1.*rand()/RAND_MAX; //lam[j]*x[j]*(1 - x[j]);
      x[j] = 1.*rand()/RAND_MAX; //lam[j]*y[j]*(1 - y[j]);
      //printf("%d %g %g\n", j, x[j], y[j]);
    }
    hs2_add(hs, x, y, 1, 1., HIST_VERBOSE);
  }
  hs2_save(hs, "HIST", wflags);

  /* now try to load histogram */
  if (0 != hs2_load(hs, "HIST", HIST_VERBOSE/*|HIST_ADDITION*/)) {
    fprintf(stderr, "cannot load histogram\n");
    return -1;
  }
  /* write again */
  hs2_save(hs, "HIST2", wflags);
  hs2_close(hs);
  return 0;
}

