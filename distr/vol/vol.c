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
  int i, t, ntot = nsteps + nequil;
  real u, k, p, s;
  real bet = 1.0/tp, fb = 0.f, udb, bvir;
  static av_t avU, avK0, avK, avp, avfb;

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
      lj_volmove(lj, volamp, tp, pressure);
    }

    if (t >= nequil) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj->rho * tp + lj->pvir);
      
      distr_add(d, lj->epot, fb, udb, 1.0);
      if ((t + 1) % nstdb == 0)
        distr_add(db, lj->epot, fb, udb, 1.0);
      av_add(&avfb, fb);
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  fb = av_getave(&avfb);
  printf("U/N = %6.3f, K/N = %6.3f, p = %6.3f, delb = %6.3f\n", u, k, p, fb);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
}

/* perform integral identity */
static void doii(distr_t *d)
{
  distr_save(d, fnds);
  distr_close(d);
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
  
  doii(d);
  doii(db);
  return 0;
}

