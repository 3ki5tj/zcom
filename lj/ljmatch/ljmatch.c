/* try to optimize the square-well potential to mimic the Lennard-Jones potential */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_AV
#define ZCOM_TMH
#include "zcom.h"

int n = 108;
int d = 3;
real rho = 0.8f;
real rcdef = 2.5f;
real tp = 1.0f;
real amp = 0.02f;
int nsteps = 100000;
int usesq = 1;
real ra = 1.0f, rb = 1.05f;
av_t avU;
double tequil = 10000;

/* tmh stuff */
double tmh_erg0, tmh_erg1, tmh_derg;
double tmh_emin, tmh_emax, tmh_de;
double tmh_ensexp = 2.0;
int tmh_dhdeorder = 1;

/* equilibrate the system */
static void equil(lj_t *lj, double tmax)
{
  double t, bet = 1.0/tp;

  for (t = 0; t < tmax; t++)
    lj_metro3d(lj, amp, bet);
}

/* match square well and Lennard-Jones */
static void match(lj_t *lj)
{
  int it;
  double t, bet = 1.0/tp;
  real esq, elj, eljn;

  elj = lj_energylj3d(lj, NULL, NULL);
  esq = lj_energysq3d(lj);
  for (t = 0; t < 10000; t++) {
    for (it = 0; it < 10; it++)
      lj_metro3d(lj, amp, bet);
    eljn = lj_energylj3d(lj, NULL, NULL);
    printf("%g %g %g %g %g\n", t, eljn - elj, lj->epot - esq, elj, esq);
    elj = eljn;
    esq = lj->epot;
  }
}

int main(void)
{
  lj_t *lj;
  tmh_t *m;

  lj = lj_open(n, d, rho, rcdef);
  lj_initsq(lj, ra, rb); /* switch to square-well */
  equil(lj, tequil);
/*
  m = tmh_open(tp, tp + 0.1, 0.01, tmh_erg0, tmh_erg1, tmh_derg,
      tmh_emin, tmh_emax, tmh_de, tmh_ensexp, tmh_dhdeorder);
*/
  match(lj);

  lj_close(lj);
  return 0;
}
