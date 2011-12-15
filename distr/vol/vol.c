/* integral identity for a volume distribution */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_DISTR
#define ZCOM_CFG
#define ZCOM_AV
#include "zcom.h"

const char *fncfg = "vol.cfg";
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

real volamp = 0.01; /* percentage */
double vmin = 400, vmax = 2000, vdel = 0.1;

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
  cfg_add(cfg, "pressure", "%r", &pressure, "pressure");
  cfg_add(cfg, "volamp", "%r", &volamp, "volume move amplitude, as a fraction");
  cfg_add(cfg, "vmin", "%lf", &vmin, "minimal volume");
  cfg_add(cfg, "vmax", "%lf", &vmax, "maximal volume");
  cfg_add(cfg, "vdel", "%lf", &vdel, "volume interval");
  cfg_add(cfg, "nstdb", "%d", &nstdb, "every # of steps to deposit to the small database (B)");
  cfg_add(cfg, "halfwin", "%d", &halfwin, "half of the number of bins in each side of the window, "
      "for the fractional identity; 0: guess, -1: adaptive");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, "
      "for mean force; 0: single bin, < 0: plain average");
  cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE);
  cfg_close(cfg);
}

/* perform a simulation, save data to ds */
static void simul(distr_t *d, distr_t *db)
{
  lj_t *lj;
  int t, ntot = nsteps + nequil, vtot = 0, vacc = 0;
  real u, k, p, vol, fv, dv, bet = 1.f/tp;
  static av_t avU, avK, avV, avp, avfv;

  lj = lj_open(N, 3, rho, rc);
  lj_initsw(lj, rs);
  printf("rc %g, rs %g, box %g\n", lj->rc, lj->rs, lj->l);
  if (initload) {
    lj_readpos(lj, lj->x, lj->v, fnpos);
    lj_force(lj);
  }

  for (t = 0; t < ntot; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    lj_vrescale(lj, tp, thermdt);
    if ((t + 1) % 10 == 0) {
      vtot += 1;
      vacc += lj_volmove(lj, volamp, tp, pressure);
    }

    if (t >= nequil) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avV, lj->vol);

      dv = lj_vir2sw3d(lj); /* r r : grad grad U */
      fv = lj->n / lj->vol - bet * lj->vir / (3*lj->vol) - bet * pressure;
      dv = -lj->n / (lj->vol * lj->vol) + bet * (2*lj->vir - dv) / (9*lj->vol*lj->vol);

      distr_add(d, lj->vol, fv, dv, 1.0);
      if ((t + 1) % nstdb == 0)
        distr_add(db, lj->vol, fv, dv, 1.0);
      av_add(&avp, lj->rho * tp - lj->vir / (3*lj->vol));
      av_add(&avfv, fv);
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  vol = av_getave(&avV)/N;
  fv = av_getave(&avfv);
  printf("U/N %6.3f, K/N %6.3f, V/N %6.3f, p %6.3f, delfv %6.3f, vacc %g%%\n",
     u, k, vol, p, fv, 100.*vacc/vtot);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
}

double LOG0 = -300;

/* compute the free energy at p + dp, also compute average volume there */
static double reweight(distr_t *d, double tp, double dp, double *vav, int raw)
{
  int i, n;
  double betdp = dp/tp, x, lnx, vol, tot = 0, lnz = LOG0, lnzv = LOG0;

  n = raw ? d->n : d->n + 1;
  for (i = 0; i < n; i++) {
    if (raw) {
      vol = d->xmin + i * d->dx;
      x = d->rho[i];
    } else {
      vol = d->xmin + (i + .5) * d->dx;
      x = d->arr[i].s;
    }
    if (x <= 0) continue;
    tot += x;
    lnx = log(x);
    lnz = lnadd(lnz, lnx - betdp * vol);
    lnzv = lnadd(lnzv, lnx - betdp * vol + log(vol));
  }
  *vav = exp(lnzv - lnz);
  return -tp * (lnz - log(tot));
}

/* do reweighting */
static void dorew(distr_t *d, distr_t *db, const char *fn)
{
  int i;
  double p, dp, vav[4], dfe[4];
  FILE *fp;

  xfopen(fp, fn, "w", return);
  for (p = 0.1; p < 0.2; p += 0.001) {
    dp = p - pressure;
    dfe[0] = reweight(db, tp, dp, &vav[0], 0);
    dfe[1] = reweight(db, tp, dp, &vav[1], 1);
    dfe[2] = reweight(d,  tp, dp, &vav[2], 0);
    dfe[3] = reweight(d,  tp, dp, &vav[3], 1);
    fprintf(fp, "%.5f ", pressure + dp);
    for (i = 0; i < 4; i++)
      fprintf(fp, "%.14e %.14e ", dfe[i], vav[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* perform integral identity */
static void doii(distr_t *d, const char *fn)
{
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

  distr_save(d, fn);
}

int main(void)
{
  distr_t *d, *db;

  loadcfg(fncfg);
  d = distr_open(vmin, vmax, vdel);
  db = distr_open(vmin, vmax, vdel);
  if (initload) { /* load previous data */
    die_if(0 != distr_load(d, fnds), "failed to load data from %s\n", fnds);
    die_if(0 != distr_load(db, fndsb), "failed to load data from %s\n", fndsb);
  }
  if (dosimul)
    simul(d, db);

  doii(d, fnds);
  doii(db, fndsb);
  dorew(d, db, "fe.dat");
  distr_close(d);
  distr_close(db);
  return 0;
}

