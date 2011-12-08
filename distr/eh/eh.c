#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_DISTR
#define ZCOM_CFG
#define ZCOM_AV
#include "zcom.h"

const char *fncfg = "eh.cfg";
const char *fnds = "ds.dat", *fndsb = "dsb.dat";
const char *fnpos = "lj.pos";
int initload = 1;

int N = 256;
real rho = 0.8f;
real rc = 3.0f;
real tp = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;

int usesw = 1; /* must use switched potential */
real rs = 2.0f;

int dosimul = 1; /* do a simulation */
int nsteps = 10000;

double emin = -1380, emax = -1240, edel = 0.2;

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
  cfg_add(cfg, "n", "%d", &N, "number of particles");
  cfg_add(cfg, "rho", "%r", &rho, "density");
  cfg_add(cfg, "rc", "%r", &rc, "cutoff distance");
  cfg_add(cfg, "rs", "%r", &rs, "switch distance");
  cfg_add(cfg, "tp", "%r", &tp, "temperature");
  cfg_add(cfg, "emin", "%lf", &emin, "minimal energy");
  cfg_add(cfg, "emax", "%lf", &emax, "maximal energy");
  cfg_add(cfg, "edel", "%lf", &edel, "energy interval");
  cfg_add(cfg, "halfwin", "%d", &halfwin, "half of the number of bins in each side of the window, for II");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, for II mean force");
  cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE);
  cfg_close(cfg);
}

/* perform a simulation, save data to ds */
static void simul(distr_t *d, distr_t *db)
{
  lj_t *lj;
  int t;
  real u, k, p;
  real bet = 1.0/tp, fb = 0.f, udb, bvir;
  static av_t avU, avK, avp, avfb;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g\n", rc, rs);
  }
  if (initload) {
    lj_readpos(lj, lj->x, lj->v, fnpos);
    lj_force(lj);
    if (usesw) fb = lj_bconfsw3d(lj, &udb, &bvir) - bet;
  }

  for (t = 0; t < nsteps; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    lj_vrescale(lj, tp, thermdt);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj->rho * tp + lj->pvir);
      if (usesw) {
        fb = lj_bconfsw3d(lj, &udb, &bvir) - bet;
        distr_add(d, lj->epot, fb, udb, 1.0);
        if ((t + 1) % 10 == 0)
          distr_add(db, lj->epot, fb, udb, 1.0);
        av_add(&avfb, fb);
      }
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

int main(void)
{
  distr_t *d, *db;

  loadcfg(fncfg);
  d = distr_open(emin, emax, edel);
  db = distr_open(emin, emax, edel);
  if (initload) { /* load previous data */
    die_if(0 != distr_load(d, fnds), "failed to load data from %s\n", fnds);
    die_if(0 != distr_load(db, fndsb), "failed to load data from %s\n", fndsb);
  }
  if (dosimul)
    simul(d, db);
  if (iitype == 0) {
    distr_aj(d, halfwin);
  } else {
    if (mfhalfwin > 0) distr_mfii(db, mfhalfwin); else distr_mf0(db);
    if (halfwin >= 0) {
      distr_ii0(db, halfwin);
    } else {
      distr_iiadp(db);
    }
  }
  distr_save(d, fnds);
  distr_close(d);
  distr_save(db, fndsb);
  distr_close(db);
  return 0;
}

