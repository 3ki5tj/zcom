#include <stdio.h>
#include "rng.h"
#include "is2.h"

#define DATAFILE "is.dat"

/* randomly pick a site and flip it */
static void mc(is_t *is, int steps, double beta)
{
  int t, id, nd;
  double p[5], acc, tot;
  double e, s, se, se2, eav, cv;
  double eref, cvref, lnzref;

  if (0 != is2_load(is, DATAFILE))
    fprintf(stderr, "cannot load prev. %s\n", DATAFILE);

  p[0] = exp(-8*beta);
  p[1] = exp(-4*beta);
  p[2] = p[3] = p[4] = 1.0;
  acc = tot = 0.0;
  s = se = se2 = 0.0;
  for (t = 1; t <= steps; t++) {
    id = is2_pick(is, &nd);
    tot += 1.0;
    if (nd >= 2 || rnd0() < p[nd]) {
      is2_flip(is, id, nd);
      acc += 1.0;
    }
    if (t % 10 == 0) {
      s += 1.0;
      se += e = is->E;
      se2 += e*e;
    }
  }
  eav = se/s;
  cv = (beta*beta)*(se2/s - eav*eav);
  lnzref = is2_exact(is, beta, &eref, &cvref);
  printf("ar: %g, eav: %.6f (%.6f), cv: %.3f (%.3f), lnz: %.6f\n", 
      acc/tot, eav, eref, cv, cvref, lnzref);
  is2_save(is, DATAFILE);
  mtsave(NULL);
}

int main(void)
{
  is_t *is;

  if ((is = is2_open(32)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  mc(is, 100000000, 0.4);
  printf("E = %d, M = %d\n", is->E, is->M);
  is2_close(is);
  return 0;
}

