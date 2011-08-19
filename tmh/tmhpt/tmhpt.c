/* average constrained
 * multi temperature */

#include <stdio.h>
#include <math.h>

#define PT2_LB 5
#define PT2_Q  10
#define L (1<<PT2_LB)
#define Q PT2_Q
#define EMIN (-2*L*L)
#define EMAX 0

#define ZCOM_PICK
#define ZCOM_RNG
#define ZCOM_TRACE
#define ZCOM_POTTS2
#define ZCOM_TMH
#include "zcom.h"

const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis";

/* regular metropolis move */
static int mcmove(potts_t *pt)
{
  int id, nb[PT2_Q], de, so, sn;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn);
  de = nb[so] - nb[sn];
  if (de <= 0 || mtrand() <= pt->uproba[de]) {
    PT2_FLIP(pt, id, so, sn, nb);
    return 1;
  } else { 
    return 0;
  }
}

/* constant temperature MC run */
static double mcrun(potts_t *pt, double tp, double *Edev,
    int tequil, int tmax)
{
  int t;
  double Esm = .0, E2sm = 0., beta = 1.0/tp;

  PT2_SETPROBA(pt, beta);
  for (t = 0; t < tequil; t++) {
    mcmove(pt);
  }
  for (t = 0; t < tmax; t++) {
    mcmove(pt);
    Esm += pt->E;
    E2sm += 1.*pt->E*pt->E;
  }
  Esm /= tmax;
  *Edev = sqrt(E2sm/tmax - Esm*Esm);
  return Esm;
}

/* energy move under a biased potential */
static int tmhmove(tmh_t *tmh, potts_t *pt)
{
  int id, so, sn, eo, en, de, nb[PT2_Q];
  double dh;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  eo = pt->E;
  en = eo + de;
  dh = tmh_hdif(tmh, en, eo); /* modified Hamiltonian */
  if (dh <= 0 || rnd0() < exp(-dh/tmh->tp)) {
    PT2_FLIP(pt, id, so, sn, nb);
    return 1;
  } else {
    return 0;
  }
}

int tmhrun(tmh_t *tmh, potts_t *pt, double trun, double t)
{
  int it;
  double amp; 
  const double ampmax = 2e-7, ampc = 0.1, lgvdt = 1e-5;

  amp = ampmax;
  tmh_settp(tmh, tmh->tp1 - 1e-6);
  for (it = 1; t <= trun; t++, it++) {
    tmhmove(tmh, pt);
    tmh_eadd(tmh, pt->E);
    tmh_dhdeupdate(tmh, pt->E, amp);

    if (it % 10 == 0) {
      tmh_tlgvmove(tmh, pt->E, lgvdt);
      /* update amplitude */
      if ((amp=ampc/t) > ampmax) amp = ampmax;
  
      if ((int)fmod(t, 100000) == 0)
        wtrace_buf("%g %d %g %g\n", t, pt->E, tmh->tp, tmh->dhde[tmh->iec]);
    }
  }
  tmh_save(tmh, fntp, fnehis, fndhde, amp, t);
  return 0;
}

#define CONFIG "pt.dat"
int main(void)
{
  potts_t *pt;
  tmh_t *tmh;
  double x, erg0, edev0, erg1, edev1, tp0 = 0.67, tp1 = 0.77;
  double emin = EMIN, emax = EMAX, de = 16;
  double amp, t0;
  int tequil = 200000, tmcrun = 2000000;
  double trun = 1000000*100;
  int initload = 0, dhdeorder = 1;

  pt = pt2_open(L, Q);
  if (initload) {
    if (pt2_load(pt, CONFIG) != 0)
      return -1;
    if (0 != tmh_loaderange(fndhde, &erg0, &erg1, &emin, &emax, &de))
      return -1;
  } else {
    /* determine the energies at the two end temperatures */
    erg0 = mcrun(pt, tp0, &edev0, tequil, tmcrun);
    printf("LowT:  tp = %g, Eav = %g (%g), Edev = %g (%g)\n", 
        tp0, erg0, erg0/pt->n, edev0, edev0/pt->n);
    erg1 = mcrun(pt, tp1, &edev1, tequil, tmcrun);
    printf("HighT: tp = %g, Eav = %g (%g), Edev = %g (%g)\n", 
        tp1, erg1, erg1/pt->n, edev1, edev1/pt->n);

    x = (erg1 - erg0)*.10;
    erg0 += x;
    erg1 -= x;  
  }

  tmh = tmh_open(tp0, tp1, erg0, erg1, EMIN, EMAX, de, 2.0, dhdeorder);

  printf("erange (%g, %g), active (%g, %g)\n", 
      tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);

  if (initload) {
    if (tmh_load(tmh, fnehis, fndhde, &amp, &t0) != 0) {
      fprintf(stderr, "cannot load tmh\n");
      return -1;
    }
    printf("continue from t = %g\n", t0);
  } else {
    t0 = 1.;
  }

  tmhrun(tmh, pt, trun, t0);

  pt2_save(pt, CONFIG);
  pt2_close(pt);
  wtrace_buf(NULL); /* finish tracing */
  tmh_close(tmh);
  mtsave(NULL);
  return 0;
}

