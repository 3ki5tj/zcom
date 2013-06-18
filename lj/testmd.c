#include "lj.h"
#include "ljmd.h"
#include "include/av.h"

const char *fnpos = "lj.pos";
int initload = 1;

int N = 108;
real rho = 0.7f;
real rc = 1000.0f; /* use half-box cut-off */
real tp = 1.5f;
real pressure = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;
int nsteps = 200000;

int usesw = 0;
real rs = 2.0f;


/* see if force matches energy */
static void foo(lj_t *lj)
{
  int i, n = lj->n;
  real f2 = 0.f, invf2, e1, e2, del = 0.01f;
  rv3_t *x = (rv3_t *) lj->x, *f = (rv3_t *) lj->f;

  for (i = 0; i < 3*n; i++) lj->x[i] += 0.01f * (2.*rnd0() - 1.);
  lj_force(lj);
  e1 = lj->epot;
  for (i = 0; i < n; i++) f2 += rv3_sqr(f[i]);
  invf2 = del / (f2 * lj->l);
  for (i = 0; i < n; i++) rv3_sinc(x[i], f[i], invf2);
  lj_force(lj);
  e2 = lj->epot;
  printf("e1 %g, e2 %g, del %g\n", e1, e2, (e1 - e2)/del);
}


int main(void)
{
  lj_t *lj;
  int t, vtot = 0, vacc = 0, isrun = 0;
  real u, k, p, rho1, bc = 0.f;
  static av_t avU, avK, avp, avbc, avrho;
  hist_t *hsvol;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g, box %g\n", rc, rs, lj->l*.5f);
  }
  if (initload && 0 != lj_readpos(lj, lj->x, lj->v, fnpos, LJ_LOADBOX)) {
    fprintf(stderr, "error loading previous coordinates from %s\n", fnpos);
  }
  foo(lj);

  hsvol = hs_open(1, 0, 5.*lj->n/lj->rho, 2.f);

  md_shiftcom(lj->v, lj->n, lj->d);
  md_shiftang(lj->x, lj->v, lj->n, lj->d);

  for (t = 1; t <= nsteps; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    lj_vrescale(lj, tp, thermdt);

    isrun = (t > nsteps / 2);
    if (t % 5 == 0) {
      /* different barostats */
/*
      vacc += lj_mctp(lj, 0.05, tp, pressure, 0, 1e300, 0, 0);
      vacc += lj_mcp(lj, 0.05, tp, pressure, 0, 1e300, 0, 0);
      lj_langtp0(lj, 1e-5, tp, pressure, 0);
      lj_langp0(lj, 1e-5, tp, pressure, 0);
*/
      lj_langp0(lj, 1e-5, tp, pressure, 0);
      vtot++;
      if (isrun) {
        av_add(&avrho, lj->n / lj->vol);
        hs_add1ez(hsvol, lj->vol, HIST_VERBOSE);
      }
    }

    if (isrun) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj_calcp(lj, tp)); /* or lj_calcpk(lj) */
      if (usesw) {
        bc = lj_bconfsw3d(lj, NULL);
        av_add(&avbc, bc);
      }
    }
    if (t % 1000 == 0) printf("t %d\n", t);
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  bc = av_getave(&avbc);
  rho1 = av_getave(&avrho);
  printf("U/N %6.3f, K/N %6.3f, p %6.3f, bc %6.3f, rho %6.3f, vacc %g%%\n",
      u, k, p, bc, rho1, 100. * vacc / vtot);
  hs_save(hsvol, "volmd.his", HIST_NOZEROES);
  hs_close(hsvol);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  return 0;
}
