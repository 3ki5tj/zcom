/* Basic molecular dynamics simulation */
#include "lj.h"
#include "ljmd.h"
#include "ljeos.h"
#include "include/av.h"



int n = 108;
real rho = 0.7f;
real rc = 1000.0f; /* use half-box cut-off */
real tp = 1.5f;
real mddt = 0.002f;
real thermdt = 0.01f;
int nequil = 10000;
int nsteps = 100000;



int main(void)
{
  lj_t *lj;
  int t;
  real u, k, p;
  double uref, pref;
  static av_t avU[1], avK[1], avp[1];

  lj = lj_open(n, 3, rho, rc);

  for (t = 1; t <= nequil + nsteps; t++) {
    lj_vv(lj, mddt);
    lj_vrescalex(lj, tp, thermdt);

    if ( t <= nequil ) continue;
    av_add(avU, lj->epot);
    av_add(avK, lj->ekin);
    av_add(avp, lj_calcp(lj, tp)); /* or lj_calcpk(lj) */
  }
  u = av_getave(avU)/n;
  k = av_getave(avK)/n;
  p = av_getave(avp);
  uref = lj_eos3dx(rho, tp, &pref, NULL, NULL, LJEOS_PVEhBHKN);
  printf("rho %g, tp %g, ep %6.3f/%g, ek %6.3f, p %6.3f/%g\n",
      rho, tp, u, uref, k, p, pref);
  lj_writepos(lj, lj->x, lj->v, "lj.pos");
  lj_close(lj);
  return 0;
}
