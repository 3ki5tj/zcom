#ifndef LB
#define LB        5
#endif
#define L         (1 << LB)
#define DATAFILE  "is.dat"

#define IS2_LB  LB
#include "ising2.h"  /* swap with the #define LB line to test two different versions */



/* randomly pick a site and flip it */
static void runmetro(ising_t *is, double steps, double beta, int ncheck)
{
  double t, acc, tot;
  double e, s1, se, se2, eav, cv;
  double eref, cvref, lnzref;
  int h, nt;
  unsigned id = 0;

  IS2_SETPROBA(is, beta);
  acc = tot = 1e-8;
  s1 = se = se2 = 0.0;
  nt = ncheck;
  for ( t = 1.0; t <= steps; t += 1.0 ) {
#ifdef IS2FUNC
    id = is2_pick(is, &h);
#else
/*
    IS2_PSEQ(is, id, h);
*/
    IS2_PICK(is, id, h);
#endif

    if (h <= 0 || mtrand() < is->uproba[h]) {
      IS2_FLIP(is, id, h);
    }
    if ( --nt == 0 ) {
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
}



/* Wolff algorithm */
static void runwolff(ising_t *is, double steps, double beta)
{
  double t, acc, tot, padd;
  double e, se, se2, eav, cv;
  double eref, cvref, lnzref;

  padd = 1 - exp(-2*beta);
  acc = tot = 1e-8;
  se = se2 = 0.0;
  for (t = 1.0; t <= steps; t += 1.0) {
    is2_wolff(is, padd);
    se += e = is->E;
    se2 += e * e;
  }
  eav = se/t;
  cv = (beta*beta)*(se2/t - eav*eav);
  lnzref = is2_exact(is, beta, &eref, &cvref);
  printf("ar: %g, eav: %.6f (%.6f), cv: %.3f (%.3f), lnz: %.6f\n",
      acc/tot, eav, eref, cv, cvref, lnzref);
}



int main(int argc, char **argv)
{
  ising_t *is;
  int method = 0, ncheck = 10;
  double tp = 2.269;

  if ((is = is2_open(L)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }

  /* try to load the previous spin configuration */
  is2_load(is, DATAFILE);

  if ( argc > 1 ) {
    method = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    tp = atof( argv[2] );
  }
  if ( argc > 0 ) {
    ncheck = atoi( argv[3] );
  }

  is = is2_open(L);
  if ( method == 0 ) {
    runmetro(is, 1e8, 1/tp, 10);
  } else if ( method == 1 ) {
    runwolff(is, 1e5, 1/tp);
  }

  printf("E = %d, M = %d\n", is->E, is->M);
  is2_em(is);
  printf("E = %d, M = %d (ref)\n", is->E, is->M);

  is2_save(is, DATAFILE);
  mtsave(NULL);
  is2_close(is);
  return 0;
}

