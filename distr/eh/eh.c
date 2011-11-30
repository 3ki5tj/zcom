#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_DISTR
#define ZCOM_AV
#include "zcom.h"

const char *fnpos = "lj.pos";
int initload = 1;

int N = 256;
real rho = 0.8f;
real rc = 3.0f;
real tp = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;

int usesw = 1;
real rs = 2.0f;

int main(void)
{
  lj_t *lj;
  int t, nsteps = 10000;
  real u, k, p;
  real bc = 0.f, udb, bvir;
  static av_t avU, avK, avp, avbc;
  distr_t *ds;

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

  ds = distr_open(-7.0*N, -3.0*N, 0.01*N);
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
  distr_save(ds, "ds.dat");
  distr_close(ds);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  return 0;
}
