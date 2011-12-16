/* integral identity for a radial distribution function */
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
real rho = 0.8f;
real rc = 3.5f;
real tp = 1.24f;
real pressure = 0.115f;
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
  cfg_add(cfg, "rmax", "%lf", &rmax, "maximal r");
  cfg_add(cfg, "rdel", "%lf", &rdel, "r interval");
  cfg_add(cfg, "halfwin", "%d", &halfwin, "half of the number of bins in each side of the window, "
      "for the fractional identity; 0: guess, -1: adaptive");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, "
      "for mean force; 0: single bin, < 0: plain average");
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
    vir1 = .5f * dfdx / dr;

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
      if (prk->in) {
        rdot = rv3_dot(prk->dx, dx); /* sign doesn't matter, to be squared */
        vir2 += (prk->psi * (rdot*rdot)/dr2 + prk->phi)/4;
      }
      kpr = getpairindex(j, k, n);
      prk = lj->pr + kpr;
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
  lj->usesw |= 0x100; /* count out-of-range pairs */
  printf("rc %g, rs %g, box %g\n", lj->rc, lj->rs, lj->l);
  if (initload) {
    lj_readpos(lj, lj->x, lj->v, fnpos);
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

/* perform integral identity */
static void doii(distr_t *d, const char *fn)
{
/*
  if (iitype == 0) {
    distr_aj(d, halfwin);
  } else if (iitype == 1 || iitype == 2) {
    if (mfhalfwin == 0) {
      distr_mf0(d);
    } else if (mfhalfwin > 0) {
      distr_mfii(d, mfhalfwin);
    } else if (mfhalfwin < 0) {
      distr_mfav(d, -mfhalfwin);
    }

    if (iitype == 1) {
      if (halfwin >= 0)
        distr_ii0(d, halfwin);
      else
        distr_iiadp(d, 1e2);
    } else {
      distr_iimf(d);
    }
  }
*/
  distr_save(d, fn);
}

int main(void)
{
  distr_t *d, *db;

  loadcfg(fncfg);
  d = distr_open(0, rmax, rdel);
  db = distr_open(0, rmax, rdel);
  if (initload) { /* load previous data */
    die_if(0 != distr_load(d, fnds), "failed to load data from %s\n", fnds);
    die_if(0 != distr_load(db, fndsb), "failed to load data from %s\n", fndsb);
  }
  if (dosimul)
    simul(d, db);

  doii(d, fnds);
  doii(db, fndsb);
  distr_close(d);
  distr_close(db);
  return 0;
}

