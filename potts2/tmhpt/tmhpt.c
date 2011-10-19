/* Tempering with modified Hamiltonian on Potts model */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L - .5)
#define EMAX (0.5)
#define EDEL 1

#define ZCOM_PICK
#define ZCOM_LOG
#define ZCOM_POTTS2
#define ZCOM_TMH
#define ZCOM_AV
#define TMH_NOCHECK  /* dangerous macro */
#include "zcom.h"

const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis", *fnpos = "pt.pos";
int bstyle = 1; /* use beta instead of T */

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
  } else return 0;
}

/* constant temperature MC run */
static double mcrun(potts_t *pt, double beta, double *dev, double *ar,
    int tequil, int tmax, const char *label)
{
  int t, acc = 0;
  double erg;
  av_t eav[1] = {{0, 0, 0}};

  PT2_SETPROBA(pt, beta);
  for (t = 0; t < tequil; t++)
    mcmove(pt);
  for (t = 0; t < tmax; t++) {
    acc += mcmove(pt);
    av_add(eav, pt->E);
  }
  erg = av_getave(eav); *dev = av_getdev(eav);
  *ar = 1.0*acc/tmax;
  printf("%6s:  T %8.4f, eav = %12.3f (%8.4f), edev = %10.4f (%6.4f), ar %5.3f%%\n", 
        label, 1/beta, erg, erg/pt->n, *dev, *dev/pt->n, *ar * 100.0);
  return erg;
}

/* energy move under a biased potential */
static void tmhmove(tmh_t *m, potts_t *pt, double beta)
{
  int id, so, sn, eo, en, de, nb[PT2_Q];
  double dh;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  eo = pt->E;
  en = eo + de;
  dh = tmh_hdif(m, en, eo); /* modified Hamiltonian */
  if (dh <= 0 || rnd0() < exp(-dh*beta))
    PT2_FLIP(pt, id, so, sn, nb);
}

static int tmhrun(tmh_t *m, potts_t *pt, double trun, double t)
{
  logfile_t *log = log_open("pt.log");
  double amp, ampmax = 2e-7, ampc = 0.02, lgvdt = 1e-5, bet;
  int it;

  tmh_settp(m, m->tp1 + bstyle*1e-8);
  bet = bstyle ? m->tp : 1/m->tp;
  for (amp = ampmax, it = 1; t <= trun; t++, it++) {
    tmhmove(m, pt, bet);
    tmh_eadd(m, pt->E);
    tmh_dhdeupdate(m, pt->E, amp);

    if (it % 10 == 0) {
      tmh_tlgvmove(m, pt->E, lgvdt);
      bet = bstyle ? m->tp : 1/m->tp;

      if ((amp = ampc/t) > ampmax) amp = ampmax; /* update amplitude */
  
      if ((int) fmod(t, 100000) == 0)
        log_printf(log, "%g %d %g %g\n", t, pt->E, m->tp, m->dhde[m->iec]);
    }
  }
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  log_close(log); /* finish tracing */
  return 0;
}

int main(void)
{
  potts_t *pt;
  tmh_t *m;
  double x, erg0, edev0, erg1, edev1, ar0, ar1;
  double op0, op1, dop;
  double tp0 = 0.67, tp1 = 0.77, dtp = 0.001;
  double beta0 = 1.50, beta1 = 1.33, dbeta = -0.001;
  double emin = EMIN, emax = EMAX, de = EDEL, derg = 32, amp, t0, ensexp = 2.0;
  int tequil = 200000, tmcrun = 2000000, trun = 1000000*200;
  int initload = 0, dhdeorder = 0;

  pt = pt2_open(L, PT2_Q);
  if (bstyle) { op0 = beta0, op1 = beta1, dop = dbeta; } 
  else { op0 = tp0, op1 = tp1, dop = dtp;}
  if (initload) {
    die_if (pt2_load(pt, fnpos) != 0, "bad %s\n", fnpos);
    if (0 != tmh_loaderange(fndhde, &op0, &op1, &dop, 
          &erg0, &erg1, &derg, &emin, &emax, &de, &ensexp, &dhdeorder))
      return -1;
  } else {
    /* determine the energies at the two end temperatures */
    if (!bstyle) { beta0 = 1.0/tp0, beta1 = 1.0/tp1; }
    erg0 = mcrun(pt, beta0, &edev0, &ar0, tequil, tmcrun, "low T");
    erg1 = mcrun(pt, beta1, &edev1, &ar1, tequil, tmcrun, "high T");
    x = (erg1 - erg0)*.10;
    erg0 += x;
    erg1 -= x;
    erg0 = dblround(erg0, derg);
  }

  m = tmh_open(op0, op1, dop, erg0, erg1, derg, EMIN, EMAX, EDEL, ensexp, dhdeorder);
  printf("erange (%g, %g), active (%g, %g)\n", m->emin, m->emax, m->erg0, m->erg1);

  if (initload) {
    die_if (tmh_load(m, fnehis, fndhde, &amp, &t0) != 0,
      "cannot load tmh from %s\n", fnehis);
    printf("continue from t = %g\n", t0);
  } else t0 = 1.;

  tmhrun(m, pt, trun, t0);

  pt2_save(pt, fnpos);
  pt2_close(pt);
  tmh_close(m);
  mtsave(NULL);
  return 0;
}

