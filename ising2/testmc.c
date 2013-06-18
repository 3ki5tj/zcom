#define LB        5
#define L         (1 << LB)
#define DATAFILE  "is.dat"

#define IS2_LB  LB
#include "ising2.h"  /* swap with the #define LB line to test two different versions */

/* randomly pick a site and flip it */
static void mc(ising_t *is, double steps, double beta, int ncheck)
{
  double t, acc, tot;
  double e, s1, se, se2, eav, cv;
  double eref, cvref, lnzref;
  int h, nt;
  unsigned id = 0;

  if (0 != is2_load(is, DATAFILE))
    fprintf(stderr, "cannot load previous %s\n", DATAFILE);
  IS2_SETPROBA(is, beta);
  acc = tot = 1e-8;
  s1 = se = se2 = 0.0;
  nt = ncheck;
  for (t = 1.0; t <= steps; t += 1.0) {
/*
    IS2_PSEQ(is, id, h);
*/
    IS2_PICK(is, id, h);

    if (h <= 0 || mtrand() < is->uproba[h]) {
      IS2_FLIP(is, id, h);
    }
    if (nt-- == 0) {
      nt = ncheck;
      s1 += 1.0;
      se += e = is->E;
      se2 += e*e;
    }
  }
  eav = se/s1;
  cv = (beta*beta)*(se2/s1 - eav*eav);
  lnzref = is2_exact(is, beta, &eref, &cvref);
  printf("ar: %g, eav: %.6f (%.6f), cv: %.3f (%.3f), lnz: %.6f\n",
      acc/tot, eav, eref, cv, cvref, lnzref);
  is2_save(is, DATAFILE);
  mtsave(NULL);
}

int main(void)
{
  ising_t *is;

  if ((is = is2_open(L)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  mc(is, 1e8, 0.4, 10);
  printf("E = %d, M = %d\n", is->E, is->M);
  is2_close(is);
  return 0;
}

