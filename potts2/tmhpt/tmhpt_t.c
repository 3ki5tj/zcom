/* Tempering with modified Hamiltonian on Potts model
 * simple version */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L)
#define EMAX 0
#define ECNT (EMAX - EMIN + 1)
#define EDEL 1

#define BLOCK 10 /* number of runs */

#define ZCOM_PICK
#define ZCOM_LOG
#define ZCOM_POTTS2
#define ZCOM_TMH
#define ZCOM_AV
#include "zcom.h"

const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis", *fnpos = "pt.pos";
double tp0 = 0.67, tp1 = 0.77, dtp = 0.001;  /* for bstyle == 0 */
double erg0 = -1760, erg1 = -832, derg = 32;
double elimit = 32e9, springk = 100.0;
double tmcrun = 1000000, trun = 1000000*100, trep = 100000;
double ampmax = 1e-6, ampc = 0.1;
double lgvdt = 1e-6;
int dhdeorder = 0;

/* energy move under modified Hamiltonian */
static int tmhmove(tmh_t *m, potts_t *pt, double beta)
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

static int ezrun(tmh_t *m, potts_t *pt, double tpinit, double tmcrun, double trun)
{
  double t;
  logfile_t *log = log_open("tmhpt.tr");

  /* constant temperature run */
  for (t = 0; t < tmcrun; t++) tmhmove(m, pt, 1.0/tpinit);
  m->elimit = elimit;
  m->springk = springk;
  tmh_initwlcvg(m, ampc, ampmax, sqrt(0.1), .99, 0, WLCVG_UPDLNFC);
  tmh_settp(m, tpinit);
  for (t = 0; t < trun; t++) {
    tmhmove(m, pt, 1.0/m->tp);
    tmh_ezmove(m, pt->E, 1.0, lgvdt);
    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g %g\n", t, pt->E, m->tp, m->dhde[m->iec], m->wl->lnf);
  }
  log_close(log);
  tmh_save(m, fntp, fnehis, fndhde, m->wl->lnf, t);
  return 0;
}

int main(void)
{
  potts_t *pt = pt2_open(L, PT2_Q);
  tmh_t *m = tmh_open(tp0, tp1, dtp, erg0, erg1, derg, EMIN, EMAX, EDEL, 2.0, dhdeorder);
  printf("erange (%g, %g), active (%g, %g)\n", m->emin, m->emax, m->erg0, m->erg1);

  ezrun(m, pt, .5 * (tp0 + tp1), tmcrun, trun);

  pt2_save(pt, fnpos);
  pt2_close(pt);
  tmh_close(m);
  mtsave(NULL);
  return 0;
}
