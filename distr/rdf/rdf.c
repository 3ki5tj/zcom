/* integral identity for a radial distribution function 
 * Copyright Cheng Zhang 2010-2012 */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_DISTR
#define ZCOM_CFG
#define ZCOM_AV
#include "zcom.h"

const char *fncfg = "rdf.cfg";
const char *fnds = "ds.dat", *fndsb = "dsb.dat";
const char *fnpos = "lj.pos";
int initload = 1;

int N = 256;
real rho = 0.7f;
real rc = 3.5f;
real tp = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;
real rs = 2.5f;

int dosimul = 1; /* do a simulation */
int nsteps = 10000, nequil = 10000;
int nstdb = 100; /* # of steps to deposit to database B */
int nstrdf = 100;

double rmax = 5.0, rdel = 0.01;

int halfwin = 50; /* half window size */
int iitype = 0; /* 0: Adib-Jarsynski; 1: modulated */
int mfhalfwin = 0; /* half window size for mean force */

double gam = 1.0;
double sampmin = 100.0;
double ajR1 = 1.0; /* used by the Adib-Jarzynski identity */
int mlimit = -1;

/* load parameters from .cfg file */
static void loadcfg(const char *fncfg)
{
  cfg_t *cfg = cfg_open(fncfg);
  if (cfg == NULL) return;
  cfg_add(cfg, "initload", "%d", &initload, "load data initially");
  cfg_add(cfg, "dosimul", "%d", &dosimul, "do simulation");
  cfg_add(cfg, "nsteps", "%d", &nsteps, "number of steps");
  cfg_add(cfg, "nequil", "%d", &nequil, "number of equilibration steps");
  cfg_add(cfg, "n", "%d", &N, "number of particles");
  cfg_add(cfg, "rho", "%r", &rho, "density");
  cfg_add(cfg, "rc", "%r", &rc, "cutoff distance");
  cfg_add(cfg, "rs", "%r", &rs, "switch distance");
  cfg_add(cfg, "tp", "%r", &tp, "temperature");
  cfg_add(cfg, "nstrdf", "%d", &nstrdf, "every # of steps to compute rdf");
  cfg_add(cfg, "nstdb", "%d", &nstdb, "every # of steps to deposit to the small database (B)");
  cfg_add(cfg, "rdel", "%lf", &rdel, "r interval");
  cfg_add(cfg, "halfwin", "%d", &halfwin, "half of the number of bins in each side of the window, "
      "for the fractional identity; 0: guess, -1: adaptive");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, "
      "for mean force; 0: single bin, < 0: plain average");
  cfg_add(cfg, "gamma", "%lf", &gam, "half window amplification factor");
  cfg_add(cfg, "sampmin", "%lf", &sampmin, "minimal number of samples to estimate sig(mf)");
  cfg_add(cfg, "ajR1", "%lf", &ajR1, "repulsion radius used by the Adib-Jarzynski identity");
  cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE);
  cfg_close(cfg);
}

/* deposit pair information */
static void deposit(distr_t *d, lj_t *lj, real bet)
{
  int i, j, k, ipr, kpr, npr = lj->npr, n = lj->n;
  real *dx, dr, dr2, rdot, dfdx, vir1, vir2 = 0.f, hbox = .5f*lj->l;
  rv3_t df, *f = (rv3_t *) lj->f;

  for (ipr = 0; ipr < npr; ipr++) {
    ljpair_t *pr = lj->pr + ipr, *prk;
    i = pr->i;
    j = pr->j;
    dx = pr->dx;
    dr2 = rv3_sqr(dx);
    dr = sqrt(dr2);
    if (dr >= hbox) continue;
    rv3_diff(df, f[i], f[j]);
    dfdx = rv3_dot(df, dx);
    vir1 = .5f * dfdx / dr; /* .5 * f . rhat */

    /* compute the second-order virial */
    if (pr->in) {
      vir2 = pr->psi * dr2 + pr->phi;
    } else {
      vir2 = 0.f;
    }
    for (k = 0; k < n; k++) {
      if (i == k || j == k) continue;
      kpr = getpairindex(i, k, n);
      prk = lj->pr + kpr;
      die_if (!(i == prk->i && k == prk->j) && !(i == prk->j && k == prk->i),
          "pair info corrupted (%d, %d) vs (%d, %d)", i, k, prk->i, prk->j);
      if (prk->in) {
        rdot = rv3_dot(prk->dx, dx); /* sign doesn't matter, to be squared */
        vir2 += (prk->psi * (rdot*rdot)/dr2 + prk->phi)/4;
      }
      kpr = getpairindex(j, k, n);
      prk = lj->pr + kpr;
      die_if (!(j == prk->i && k == prk->j) && !(j == prk->j && k == prk->i),
          "pair info corrupted (%d, %d) vs (%d, %d)", j, k, prk->i, prk->j);
      if (prk->in) {
        rdot = rv3_dot(prk->dx, dx);
        vir2 += (prk->psi * (rdot*rdot)/dr2 + prk->phi)/4;
      }
    }
    distr_add(d, dr, bet*vir1, -bet*vir2, 1.0);
  }
}

/* perform a simulation, save data to ds */
static void simul(distr_t *d, distr_t *db)
{
  lj_t *lj;
  int t, ntot = nsteps + nequil, nb = 0;
  real u, k, bet = 1.f/tp;
  static av_t avU, avK;

  lj = lj_open(N, 3, rho, rc);
  lj_initsw(lj, rs);
  lj->usesw |= LJ_SWALLPAIRS; /* count out-of-range pairs */
  printf("rc %g, rs %g, box %g\n", lj->rc, lj->rs, lj->l);
  if (initload) {
    lj_readpos(lj, lj->x, lj->v, fnpos, 0);
    lj_force(lj);
  }

  for (t = 0; t < ntot; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    lj_vrescale(lj, tp, thermdt);

    if (t >= nequil) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);

      if ((t + 1) % nstrdf == 0) {
        deposit(d, lj, bet);
        if (++nb % nstdb == 0)
          deposit(db, lj, bet);
      }
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  printf("U/N %6.3f, K/N %6.3f\n", u, k);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
}

/* make an array of normalization factor */
static void mkrdfnorm(distr_t *d, double vol)
{
  double dr = d->dx, r;
  int i, max = d->n + 100;

  vol *= M_PI/6;
  xnew(d->norm, max);
  for (i = 0; i < max; i++) {
    r = i * dr;
    d->norm[i] = 4*M_PI*r*r/vol;
  }
}

/* the Adib-Jazynski identity from the original paper
 * for rdf only */
INLINE void distr_ajrdf(distr_t *d, double R1)
{
  int i1, i, j, n = d->n;
  double hbox, vol, bet = 1.0/tp, his;
  double tot, s1, den, den0, dx = d->dx;
  double u, *ua, *ub, *xp, *fr;
  real r, fscal, psi, xi, vr;
  distrsum_t *ds;
  lj_t *lj = lj_open(N, 3, rho, rc);
  lj_initsw(lj, rs);

  /* compute the total visits */
  for (tot = 0., i = 0; i < n; i++) tot += d->arr[i].s;

  i1 = (int)(R1/dx + .1);
  R1 = i1 * dx;
  vol = N/rho;
  hbox = pow(vol, 1.0/3) * .5;
  den0 = tot * (1 - R1*R1*R1/(hbox*hbox*hbox));

  /* preparation */
  xnew(xp, n + 1);
  xnew(ua, n + 1);
  xnew(ub, n + 1);
  xnew(fr, n + 1);
  for (his = 0, i = 0; i <= n; i++) {
    r = (real) ((i + .5) * d->dx);
    vr = lj_potsw(lj, r, &fscal, &psi, &xi);
    xp[i] = ua[i] = ub[i] = fr[i] = 0.;
    vr *= bet;
    if (vr > 100) vr = 100;
    fr[i] = fscal * r;
    xp[i] = exp(vr);
    d->lnrho[i] = fr[i];

    /* sum exp(bet*vr) within the volume omega* to his */
    if (i >= i1 && i < n)
      his += d->arr[i].s * xp[i];
    
    /* compute vector field exp(bet*vr)/(3r^2) * (r^3 - R1^3) */
    if (i > i1)
      ua[i] = xp[i]/(3*r*r)*(r*r*r - R1*R1*R1);

    /* compute vector field exp(bet*vr)/(3r^2) * (R1^3 - R^3) */
    ub[i] = xp[i]/(3*r*r)*(R1*R1*R1 - hbox*hbox*hbox); 
  }
  printf("his = %g den0 %g\n", his, den0);

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute u'*f */
    for (s1 = 0, j = 0; j < n; j++) {
      ds = d->arr + j;
      if (ds->s <= 0.) continue;
      if (j >= i) u = ub[j]; else u = ua[j];
      s1 += u * (ds->sf - bet * ds->s * fr[j]);
    }
    den = den0 * xp[i];
    d->his[i] = his/den;
    d->rho[i] = (his + s1)/den;
  }
  lj_close(lj);
  free(xp);
  free(ua);
  free(ub);
}

static void doii(distr_t *d, const char *fn)
{
  if (iitype == 9) {
    distr_ajrdf(d, ajR1);
  } else {
    distr_iiez(d, iitype, halfwin, mfhalfwin, gam, mlimit, sampmin);
  }
  distr_save(d, fn);
}

int main(void)
{
  distr_t *d, *db;

  loadcfg(fncfg);
  rmax = (real) ((int) (.5 * pow(N/rho, 1.0/3)/rdel) * rdel);
  printf("rmax = %g\n", rmax);
  d = distr_open(0, rmax, rdel);
  db = distr_open(0, rmax, rdel);
  if (initload) { /* load previous data */
    die_if(0 != distr_load(d, fnds), "failed to load data from %s\n", fnds);
    die_if(0 != distr_load(db, fndsb), "failed to load data from %s\n", fndsb);
  }
  mkrdfnorm(d, N/rho);
  mkrdfnorm(db, N/rho);
  if (dosimul)
    simul(d, db);

  doii(d, fnds);
  doii(db, fndsb);
  distr_close(d);
  distr_close(db);
  return 0;
}

