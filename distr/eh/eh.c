#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_DISTR
#define ZCOM_CFG
#define ZCOM_AV
#include "zcom.h"

const char *fncfg = "eh.cfg";
const char *fnds = "ds.dat";
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

double emin = -7.0, emax = -3.0, edel = 0.001;

int halfwin = 50; /* half window size */
int iitype = 0; /* 0: Adib-Jarsynski, 1: modulated */

/* load parameters from .cfg file */
static void loadcfg(const char *fncfg)
{
  cfg_t *cfg = cfg_open(fncfg);
  if (cfg == NULL) return;
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
  cfg_add(cfg, "halfwin", "%d", &halfwin, "II: number of bin in each side of the window");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE);
  cfg_close(cfg);
}

/* perform a simulation, save data to ds */
static void simul(distr_t *ds)
{
  lj_t *lj;
  int t;
  real u, k, p;
  real bc = 0.f, udb, bvir;
  static av_t avU, avK, avp, avbc;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g\n", rc, rs);
  }
  if (initload) {
    lj_readpos(lj, lj->x, lj->v, fnpos);
    lj_force(lj);
    if (usesw) bc = lj_bconfsw3d(lj, &udb, &bvir);
  }

  for (t = 0; t < nsteps; t++) {
    lj_vv(lj, mddt);
    lj_vrescale(lj, tp, thermdt);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj->rho * tp + lj->pvir);
      if (usesw) {
        bc = lj_bconfsw3d(lj, &udb, &bvir);
        distr_add(ds, lj->epot, bc, 1.0/tp, udb, 1.0);
        av_add(&avbc, bc);
      }
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  bc = av_getave(&avbc);
  printf("U/N = %6.3f, K/N = %6.3f, p = %6.3f, bc = %6.3f\n", u, k, p, bc);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
}

int main(void)
{
  distr_t *ds;

  loadcfg(fncfg);
  ds = distr_open(emin*N, emax*N, edel*N);
  if (dosimul) {
    simul(ds);
  } else { /* load previous data */
    die_if(0 != distr_load(ds, fnds), "failed to load data from %s\n", fnds);
  }
  if (iitype == 0) {
    distr_aj(ds, halfwin);
  } else {
    distr_ii0(ds, halfwin);
  }
  distr_save(ds, fnds);
  distr_close(ds);
  return 0;
}

