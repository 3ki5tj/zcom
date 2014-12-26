/* basic Monte Carlo simulation */
#include "lj.c"
#include "include/av.h"

int n = 108;
int d = 3;
real rho = (real) 0.5;
real rcdef = 2.5f;
real tp = (real) 2.0;
real pressure = (real) 1.0;
real amp = (real) 0.3;
int nsteps = 1000000;
av_t avU, avp, avrho;

int main(void)
{
  int t, acc = 0;
  lj_t *lj;
  real u, p, p1;

  lj = lj_open(n, d, rho, rcdef);
  lj->dof = n * d; /* set the DOF for MC */
  for (t = 1; t <= nsteps; t++) {
    acc += lj_metro(lj, amp, 1.0f/tp);
    if ( t < nsteps / 2 ) continue;
    av_add(&avU, lj->epot);
    av_add(&avp, lj_calcp(lj, tp));
  }
  u = av_getave(&avU)/lj->n;
  p = av_getave(&avp);
  lj_gettail(lj, rho, n, &p1);
  printf("t %d, l %g, erg %g, p %g/%g, acc %.2f%%, dof/d %g\n",
      nsteps, lj->l, u, p, p1, 100.*acc/nsteps, 1.*lj->dof/lj->d);
  lj_close(lj);
  return 0;
}

