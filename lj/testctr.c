#include "lj.c"
#include "include/av.h"

int n = 108;
int d = 3;
real rho = 0.8f;
real rc = 2.5f;
real tp = 1.0f;
real amp = 0.04f;
int nsteps = 100000;
int usesw = 1;
real rs = 1.5f; /* for switched potential */

/* local move */
static real tcr3d_l(lj_t *lj, int tmax, real amp, real umax, real udel)
{
  int i, t;
  real du, du2, betp, vir, rmin, xi[3];
  hist_t *hs = hs_open(1, -umax, umax, udel);
  av_t avdu[1];

  amp /= lj->l;
  av_clear(avdu);
  for (t = 0; t < tmax; t++) {
    i = lj_randmv3d(lj, xi, amp*sqrt(lj->n));
    du = lj_depot3d(lj, i, xi, &vir, &rmin);
    hs_add(hs, &du, 1.0, 0);
    av_add(avdu, du); 
  }
  du = av_getave(avdu);
  du2 = av_getvar(avdu);
  betp = (du2 > 0.) ? 2.f*du/du2 : 0.;
  printf("local perturbation: du %g, du2 %g, betp %g, tp %g\n", du, du2, betp, 1.0/betp);
  hs_save(hs, "de_l.his", HIST_ADDAHALF);
  return betp;
}

static real tcr3d_g(lj_t *lj, int tmax, real amp, real umax, real udel)
{
  int i, t, d, n = lj->n;
  real du, ep = lj->epot, du2, betp;
  rv3_t *x0, *x = (rv3_t *) lj->x;
  hist_t *hs = hs_open(1, -umax, umax, udel);
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
    du = lj_energy3d(lj) - ep;

    hs_add(hs, &du, 1.0, 0);
    av_add(avdu, du); 
  }
  du = av_getave(avdu);
  du2 = av_getvar(avdu);
  betp = (du2 > 0.) ? 2.f*du/du2 : 0.;
  printf("global perturbation: du %g, du2 %g, betp %g, tp %g\n", du, du2, betp, 1.0/betp);
  hs_save(hs, "de.his", HIST_ADDAHALF);
  free(x0);
  return betp;
}

int main(void)
{
  int t, acc = 0;
  lj_t *lj = lj_open(n, d, rho, rc);
  real u, bc = 1.0f/tp, bc0 = bc;

  if (usesw) lj_initsw(lj, rs);
  for (t = 0; t < nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
  }
  u = lj->epot;
  if (usesw) {
    lj_force(lj);
    bc = lj_bconfsw3d(lj, NULL);
    bc0 = lj->lap/lj->f2;
  }
  printf("finish equil., u %g, bet %g, bc %g, bc0 %g, acc %g\n",
      u/lj->n, 1.0/tp, bc, bc0, 1.*acc/nsteps);
  tcr3d_l(lj, 10000000, 0.001, 1.0, 0.01);
  tcr3d_g(lj, 1000000, 0.005, 10.0, 0.1);
  lj_close(lj);
  return 0;
}
