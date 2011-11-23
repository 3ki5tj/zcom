/* Tempering with modified Hamiltonian on Potts model */
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
int bstyle = 0; /* 1: use beta, 0: use T */
double tp0 = 0.67, tp1 = 0.77, dtp = 0.001;  /* for bstyle == 0 */
double beta0 = 1.50, beta1 = 1.33, dbeta = -0.002; /* for bstyle == 1 */
double erg0 = -1760, erg1 = -832, derg = 32, elimit = 32;
int tequil = 200000, tmcrun = 2000000;
double trun = 1000000*300, trep = 100000;
double ampmax = 1e-5, ampc = 10.0;
double lgvdt = 1e-6;
int dhdeorder = 1;
int entropic = 0; /* entropic sampling */
double entampmax = 5e-4, entampc = 5000.0;
int initload = 0; /* continue from a previous run */
int update = 1; /* update dhde */
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
    if (fabs(epot - m->ec) < elimit || epot < m->erg0 || epot > m->erg1) {
      tmh_updhde(m, (epot - m->ec) * amp);
    }
    tmh_lgvmove(m, epot, lgvdt);
    bet = bstyle ? m->tp : 1/m->tp;

    if ((int) fmod(t, trep) == 0)
      log_printf(log, "%g %d %g %g %g\n", t, pt->E, m->tp, m->dhde[m->iec], amp);
  }
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

/* compile mh according to dhde, similar to tmh_caclmh()
 * for dhdeorder == 0, assume EDEL = 1.0 */
INLINE void tmhcalcmh0(tmh_t *m)
{
#if 0
  /* simple implementation */
  int i, en = m->en;
  for (m->mh[0] = EMIN, i = 0; i < en; i++)
    m->mh[i+1] = m->mh[i] + tmh_getdhde(m, EMIN + i * EDEL);
#endif
  /* slightly faster implementation */
  int i, j, k, en = m->en, i0 = (int)(m->erg0 - EMIN + .5), di = (int)(m->derg + .5);
  double x = EMIN, dh = m->dhde[0];

  for (m->mh[i = 0] = x; i < i0; ) m->mh[++i] = (x += dh);
  for (j = 0; j < m->ergn; j++)
    for (dh = m->dhde[j], k = 0; k < di; k++)
      m->mh[++i] = (x += dh);
  for (; i < en; ) m->mh[++i] = (x += dh);
}

typedef struct {
  double ec;
  int iec;
} e2tdata_t;

e2tdata_t e2tarr[ECNT];

static void tmhcalce2tarr(tmh_t *m)
{
  int e, ie;
  for (e = EMIN; e <= EMAX; e++) {
    tmh_setec(m, e);
    ie = e - EMIN;
    e2tarr[ie].ec = m->ec;
    e2tarr[ie].iec = m->iec;
  }
}

/* static energy move */
INLINE int tmhmove0(tmh_t *m, potts_t *pt, double beta)
{
  int id, so, sn, de, nb[PT2_Q], acc;
  double dh;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  if (de > 0) { /* avoid calling hdif() if possible */
    int eo = pt->E;
    int en = eo + de;
    dh = m->mh[en - EMIN] - m->mh[eo - EMIN]; /* modified Hamiltonian */
    acc = (dh <= 0) || (rnd0() < exp(-dh*beta));
  } else acc = 1;
  if (acc) { PT2_FLIP(pt, id, so, sn, nb); return de; }
  else return 0;
}

/* entropic sampling: simple implementation */
static int tmhrun_ent0(tmh_t *m, potts_t *pt, double bet, double trun, double t,
    logfile_t *log)
{
  double amp = entampmax;
  int it, ie, de, wt = 1, ierg0 = (int)(erg0 - .5), ierg1 = (int)(erg1 - .5);

  tmh_setec(m, pt->E);

  /* production run */
  for (; t < trun; t += BLOCK) {
    amp = dblmin(entampc/(t + .1), entampmax);
    for (it = 0; it < BLOCK; it++) {
      de = tmhmove(m, pt, bet);
      if (de != 0 || pt->E < ierg0 || pt->E > ierg1) {
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

/* entropic sampling  */
static int tmhrun_ent1(tmh_t *m, potts_t *pt, double bet, double trun, double t,
    logfile_t *log)
{
  double amp = entampmax;
  int it, ie, de, wt = 1, ierg0 = (int)(erg0 - .5), ierg1 = (int)(erg1 - .5);

/* map energy to temperature, set ec, iec, tp, itp */
#define E2T(ie) { e2tdata_t *e2t = e2tarr + (ie); \
  m->ec = e2t->ec; m->iec = e2t->iec; }

  tmhcalce2tarr(m);
  E2T(pt->E - EMIN);
  tmhcalcmh0(m);

  /* production run */
  for (; t < trun; t += BLOCK) {
    amp = dblmin(entampc/(t + .1), entampmax);
    for (it = 0; it < BLOCK; it++) {
      de = tmhmove0(m, pt, bet);
      if (de != 0 || pt->E < ierg0 || pt->E > ierg1) {
        ie = pt->E - EMIN;
        m->tpehis[ie] += wt; wt = 1;
        tmh_updhde(m, (pt->E - m->ec) * amp);
        E2T(ie);
      } else wt++;
    }

    if ((int) fmod(t, 1000) == 0) {
      tmhcalcmh0(m);
      if ((int) fmod(t, trep) == 0)
        log_printf(log, "%g %d %g %g %g\n", t, pt->E, 1.0/bet, m->dhde[m->iec], amp);
    }
  }
  tmh_setec(m, pt->E);
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

unsigned mhproba[ECNT][5];

/* compile transition probability under modified transitions */
static void tmhcalcproba(tmh_t *m, double bet)
{
  int e, ie, ie1, de;
  double x;

  for (e = EMIN; e <= EMAX; e++) {
    ie = e - EMIN;
    tmh_setec(m, e);
    mhproba[ie][0] = 0xffffffffu;
    for (de = 1; de <= 4; de++) {
      ie1 = ie + de;
      if (ie1 >= ECNT) {
        mhproba[ie][de] = 0;
      } else if (m->mh[ie1] <= m->mh[ie]) {
        mhproba[ie][de] = 0xffffffffu;
      } else {
        x = exp(-bet * (m->mh[ie1] - m->mh[ie]));
        mhproba[ie][de] = (unsigned)(0xffffffffu * x);
      }
    }
  }
}

/* static energy move, for entropic sampling only
 * assumes a static dhde */
INLINE int tmhmove1(potts_t *pt)
{
  int id, so, sn, de, nb[PT2_Q];

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  if (de <= 0 || mtrand() <= mhproba[pt->E - EMIN][de]) {
    PT2_FLIP(pt, id, so, sn, nb); return de;
  } else return 0;
}

/* entropic sampling: constant dhde, only for verify */
static int tmhrun_entvrf(tmh_t *m, potts_t *pt, double bet, double trun, double t,
    logfile_t *log)
{
  int it;

  tmhcalcmh0(m);
  tmhcalcproba(m, bet);
  /* production run */
  for (; t < trun; t += BLOCK) {
    for (it = 0; it < BLOCK; it++) tmhmove1(pt);
    if ((int) fmod(t, trep) == 0) log_printf(log, "%g %d\n", t, pt->E);
  }
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
    tmh_savedhde(m, "a.e", amp, t0);
    opinit = m->tp;
    printf("continue from t %g, op %g\n", t0, opinit);
  } else t0 = 0.;

  if (entropic) {
    die_if (dhdeorder != 0, "must use zeroth order for entropic sampling\n");
    if (update) {
      //tmhrun_ent1(m, pt, beta0, trun, t0, log);
      tmhrun_ent0(m, pt, beta0, trun, t0, log);
    } else {
      tmh_loaddhde(m, fndhde, &amp, &t0);
      tmhrun_entvrf(m, pt, beta0, trun, t0, log);
    }
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

