#include "lj.c"
#include "include/av.h"
#include "include/hist.c"

int n = 108;
int d = 3;
real rho = 0.8f;
real rc = 2.5f;
real tp = 1.0f;
real amp = 0.04f;
int nsteps = 10000;
int usesq = 0, usesw = 1;
real ra = 1.0f, rb = 1.25f; /* for square well potential */
real rs = 1.5f; /* for switched potential */

static real tcr3d(lj_t *lj, int tmax, real amp, real umax, real du)
{
  int i, t, d, n = lj->n;
  real u, ep = lj->epot, dv, bet;
  rv3_t *x0, *x = (rv3_t *) lj->x;
  hist_t *hs = hs_open(1, -umax, umax, du);
  av_t avdu[1];

  xnew(x0, n);
  rv3_ncopy(x0, x, n);
  amp /= lj->l;
  av_clear(avdu);
  for (t = 0; t < tmax; t++) {
    /* randomly displace x0 */
    for (i = 0; i < n; i++)
      for (d = 0; d < 3; d++)
        x[i][d] = x0[i][d] + amp * (2 * rnd0() - 1);
    u = lj_energy3d(lj);
    du = u - ep;
    hs_add(hs, &du, 1.0, 0);
    av_add(avdu, du); 
  }
  du = av_getave(avdu);
  dv = av_getdev(avdu);
  bet = (dv > 0.) ? 2.f*du/(dv*dv) : 0.;
  printf("du %g, dev %g, bet %g, tp %g\n", du, dv, bet, 1.0/bet);
  hs_save(hs, "de.his", HIST_ADDAHALF);
  free(x0);
  return bet;
}

int main(void)
{
  int t, acc = 0;
  lj_t *lj = lj_open(n, d, rho, rc);
  real u, bc = 1.0f/tp;

  if (usesq) lj_initsq(lj, ra, rb);
  else if (usesw) lj_initsw(lj, rs);
  for (t = 0; t < nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
  }
  u = lj->epot;
  if (usesw) {
    lj_force(lj);
    bc = lj_bconfsw3d(lj, NULL);
  }
  printf("finish equil., u %g, bet %g, bc %g, acc %g\n",
      u/lj->n, 1.0/tp, bc, 1.*acc/nsteps);
  tcr3d(lj, 100, 0.01, 20.0, 0.1);
  lj_close(lj);
  return 0;
}
