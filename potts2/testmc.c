#include <stdio.h>
#include "rng.c"
#include "util.h"

#define LB        5
#define L         (1 << LB)
#define Q         10
#define DATAFILE  "pt.dat"

#define PT2_LB  LB 
#define PT2_Q   Q
#include "potts2.c"  /* swap with the #define LB line to test two different versions */

/* randomly pick a site and flip it */
static void mc(potts_t *pt, double steps, double beta, int ncheck)
{
  double t, acc, tot;
  double e, s1, se, se2, eav, cv;
  int nb[Q], de, nt, so, sn;
  unsigned id = 0;

  if (0 != pt2_load(pt, DATAFILE))
    fprintf(stderr, "cannot load previous %s\n", DATAFILE);
  PT2_SETPROBA(pt, beta);
  acc = tot = 1e-8;
  s1 = se = se2 = 0.0;
  nt = ncheck;
  for (t = 1.0; t <= steps; t += 1.0) {
    /* regular metropolis */
/*
    PT2_PSEQ(pt, id, nb);
    //PT2_PICK(pt, id, nb);
    PT2_NEWFACE(pt, id, so, sn);
    de = nb[so] - nb[sn];
    if (de <= 0 || mtrand() <= pt->uproba[de]) {
      PT2_FLIP(pt, id, so, sn, nb);
    }
*/

    /* heat bath */

    //PT2_PSEQ(pt, id, nb);
    PT2_PICK(pt, id, nb);
    PT2_HEATBATH(pt, id, so, sn, nb);
    if (so != sn) {
      PT2_FLIP(pt, id, so, sn, nb);
    }


    if (nt-- == 0) {
      nt = ncheck;
      s1 += 1.0;
      se += e = pt->E;
      se2 += e*e;
    }
  }
  eav = se/s1;
  cv = (beta*beta)*(se2/s1 - eav*eav);
  printf("ar: %g, eav: %.6f, cv: %.3f\n", 
      acc/tot, eav, cv);
  pt2_save(pt, DATAFILE);
  mtsave(NULL);
}

int main(void)
{
  potts_t *pt;
  if ((pt = pt2_open(L, Q)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  mc(pt, 1e8, 0.5, 10);
  printf("E = %d, M = %d\n", pt->E, pt->M[0]);
  pt2_close(pt);
  return 0;
}

