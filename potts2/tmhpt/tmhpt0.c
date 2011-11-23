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
#define TMH_NOCHECK  /* dangerous macro */
#include "zcom.h"

const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis", *fnpos = "pt.pos";
int bstyle = 0; /* use beta instead of T */
double tp0 = 0.67, tp1 = 0.77, dtp = 0.001;  /* for bstyle == 0 */
double beta0 = 1.50, beta1 = 1.33, dbeta = -0.002; /* for bstyle == 1 */
double erg0 = -1760, erg1 = -832, derg = 32, elimit = 32;
int tequil = 200000, tmcrun = 2000000;
double trun = 1000000*100, trep = 100000;
double ampmax = 1e-5, ampc = 10.0;
double lgvdt = 1e-6;
int dhdeorder = 0;
int entropic = 0; /* entropic sampling */
double entampmax = 5e-4, entampc = 5000.0;
int initload = 0; /* continue from a previous run */
int guesse = 0; /* guess E range */

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

static int tmhrun(tmh_t *m, potts_t *pt, double opinit, double trun, double t,
    logfile_t *log)
{
  double bet, epot, amp;
  int it;

  tmh_settp(m, opinit);
  bet = bstyle ? m->tp : 1/m->tp;
  for (; t < trun; t += BLOCK) {
    amp = dblmin(ampc/(t + .1), ampmax);
    for (it = 0; it < BLOCK; it++) {
      tmhmove(m, pt, bet);
      tmh_eadd(m, pt->E + .5);
    }
    epot = pt->E;
    if (fabs(epot - m->ec) < elimit || epot < m->erg0 || epot > m->erg1)
      tmh_updhde(m, (epot - m->ec) * amp);
    tmh_lgvmove(m, epot, lgvdt);
    bet = bstyle ? m->tp : 1/m->tp;

    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g %g\n", t, pt->E, m->tp, m->dhde[m->iec], amp);
  }
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

/* entropic sampling */
static int tmhrun_ent0(tmh_t *m, potts_t *pt, double bet,
    double trun, double t, logfile_t *log)
{
  double amp = entampmax;
  int it, ie, de, wt = 1;

  /* production run */
  tmh_setec(m, pt->E);
  for (; t < trun; t += BLOCK) {
    amp = dblmin(entampc/(t + .1), entampmax);

    for (it = 0; it < BLOCK; it++) {
      de = tmhmove(m, pt, bet);
      if (de != 0) {
        ie = pt->E - EMIN;
        m->tpehis[ie] += wt; wt = 1;
        tmh_updhde(m, (pt->E - m->ec) * amp);
        tmh_setec(m, pt->E);
      } else wt++;
    }

    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g %g\n", t, pt->E, 1.0/bet, m->dhde[m->iec], amp);
  }
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

int main(void)
{
  potts_t *pt;
  tmh_t *m;
  double x, edev0, edev1, ar0, ar1;
  double op0, op1, dop, opinit;
  double emin = EMIN, emax = EMAX, de = EDEL, amp, t0, ensexp = 2.0;
  logfile_t *log = log_open("tmhpt.tr");

  pt = pt2_open(L, PT2_Q);
  if (bstyle) { op0 = beta0, op1 = beta1, dop = dbeta; }
  else { op0 = tp0, op1 = tp1, dop = dtp;}
  opinit = .5 * (op0 + op1);
  if (initload) {
    die_if (pt2_load(pt, fnpos) != 0, "bad %s\n", fnpos);
    if (0 != tmh_loaderange(fndhde, &op0, &op1, &dop,
          &erg0, &erg1, &derg, &emin, &emax, &de, &ensexp, &dhdeorder))
      return -1;
  } else if (guesse) {
    /* determine the energies at the two end temperatures */
    if (!bstyle) beta0 = 1.0/tp0, beta1 = 1.0/tp1;
    opinit = op1 + (bstyle ? 1e-8 : -1e-8);
    erg0 = mcrun(pt, beta0, &edev0, &ar0, tequil, tmcrun, "low T");
    erg1 = mcrun(pt, beta1, &edev1, &ar1, tequil, tmcrun, "high T");
    x = (erg1 - erg0)*.10;
    erg0 += x;
    erg1 -= x;
    erg0 = dblround(erg0, derg);
  } else {
    double betainit = bstyle ? opinit : 1.0/opinit;
    mcrun(pt, betainit, &edev0, &ar0, tequil, tmcrun, "equilibration");
  }

  m = tmh_open(op0, op1, dop, erg0, erg1, derg, emin, emax, de, ensexp, dhdeorder);
  printf("erange (%g, %g), active (%g, %g)\n", m->emin, m->emax, m->erg0, m->erg1);

  if (initload) {
    die_if (tmh_load(m, fnehis, fndhde, &amp, &t0) != 0,
      "cannot load tmh from %s\n", fnehis);
    opinit = m->tp;
    printf("continue from t %g, op %g\n", t0, opinit);
  } else t0 = 0.;

  if (entropic) {
    tmhrun_ent0(m, pt, beta0, trun, t0, log);
  } else {
    tmhrun(m, pt, opinit, trun, t0, log);
  }

  pt2_save(pt, fnpos);
  pt2_close(pt);
  tmh_close(m);
  mtsave(NULL);
  log_close(log); /* finish tracing */
  return 0;
}

