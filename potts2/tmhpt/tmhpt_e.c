/* algorihtm E on Potts model, zeroth and first order */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L)
#define EMAX 0
#define EDEL 1

#define ZCOM_PICK
#define ZCOM_LOG
#define ZCOM_POTTS2
#define ZCOM_TMH
#define TMH_NOCHECK  /* dangerous macro */
#include "zcom.h"

int dhdeorder = 1, easy = 1;
double beta0 = 1.40, beta1 = 1.33, dbeta = -0.01; /* only beta0 is actually used */
double erg0 = -1760, erg1 = -832, derg = 32;
double trun = 1000000*100, trep = 100000, tmcrun = 2000000;
double entampmax = 1e-3, entampc = 5000.0;
const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis", *fnpos = "pt.pos";

/* energy move under modified Hamiltonian */
static int move(tmh_t *m, potts_t *pt, double beta)
{
  int id, so, sn, eo, en, de, nb[PT2_Q], acc;
  double dh;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  if (de > 0) { /* avoid calling hdif() if possible */
    eo = pt->E;
    en = eo + de;
    dh = tmh_hdif(m, en, eo); /* modified Hamiltonian */
    acc = (dh <= 0) || (rnd0() < exp(-dh*beta));
  } else acc = 1;
  if (acc) { PT2_FLIP(pt, id, so, sn, nb); return de; }
  else return 0;
}

/* entropic sampling */
static int run(tmh_t *m, potts_t *pt, double tmcrun, double trun)
{
  double t, amp;
  int it, ie, de, wt = 1;
  logfile_t *log = log_open("tmhpt.tr");

  for (t = 1; t <= tmcrun; t++) /* equilibration */
    move(m, pt, beta0);

#define BLOCK 10 /* block to avoid overhead */
  tmh_setec(m, pt->E);
  for (t = 0; t < trun; t += BLOCK) { /* production */
    amp = dblmin(entampc/(t + .1), entampmax);
    
    for (it = 0; it < BLOCK; it++) {
      de = move(m, pt, beta0);
      if (de != 0) {
        ie = pt->E - EMIN;
        m->tpehis[ie] += wt; wt = 1;
        m->dhde[m->iec] = dblmax(m->dhde[m->iec] + (pt->E - m->ec) * amp, 0);
        tmh_setec(m, pt->E);
      } else wt++;
    }

    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g\n", t, pt->E, m->dhde[m->iec], amp);
  }
  log_close(log);
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

/* same as run(), but use shortcut routines */
static int ezrun(tmh_t *m, potts_t *pt, double tmcrun, double trun)
{
  double t;
  logfile_t *log = log_open("tmhpt.tr");

  for (t = 1; t <= tmcrun; t++) move(m, pt, beta0);
  tmh_initwlcvg(m, entampc, entampmax, sqrt(0.1), .95, 0, WLCVG_UPDLNFC);
  tmh_setec(m, pt->E);
  for (t = 0; t < trun; t++) { /* production */
    move(m, pt, beta0);
    tmh_ezmoves(m, pt->E, 1.0);
    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g\n", t, pt->E, m->dhde[m->iec], m->wl->lnf);
  }
  log_close(log);
  tmh_save(m, fntp, fnehis, fndhde, m->wl->lnf, t);
  return 0;
}

int main(void)
{
  potts_t *pt = pt2_open(L, PT2_Q);
  tmh_t *m = tmh_open(beta0, beta1, dbeta, erg0, erg1, derg, EMIN, EMAX + EDEL, EDEL, 2.0, dhdeorder);
  
  if (easy) ezrun(m, pt, tmcrun, trun);
  else run(m, pt, tmcrun, trun);

  pt2_save(pt, fnpos);
  pt2_close(pt);
  tmh_close(m);
  mtsave(NULL);
  return 0;
}

