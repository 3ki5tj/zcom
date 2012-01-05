#include "lj.c"
#include "include/av.h"

const char *fnpos = "lj.pos";
int initload = 1;

int N = 256;
real rho = 0.7f;
real rc = 3.0f;
real tp = 1.2f;
real mddt = 0.002f;
real thermdt = 0.01f;
real pressure = 0.1f;
int nsteps = 10000;

int usesw = 1;
real rs = 2.0f;

/* see if force matches energy */
static void foo(lj_t *lj)
{
  int i, n = lj->n;
  real f2 = 0.f, invf2, e1, e2, del = 0.1f;
  rv3_t *x = (rv3_t *) lj->x, *f = (rv3_t *) lj->f;

  for (i = 0; i < 3*n; i++) lj->x[i] += 0.01f * (2.*rnd0() - 1.);
  lj_force(lj);
  e1 = lj->epot;
  for (i = 0; i < n; i++) f2 += rv3_sqr(f[i]);
  invf2 = del/f2/lj->l;
  for (i = 0; i < n; i++) rv3_sinc(x[i], f[i], invf2);
  lj_force(lj);
  e2 = lj->epot;
  printf("e1 %g, e2 %g, del %g\n", e1, e2, (e1 - e2)/del);
}

int main(void)
{
  lj_t *lj;
  int t, vtot = 0, vacc = 0;
  real u, k, p;
  real bc = 0.f;
  static av_t avU, avK, avp, avbc;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g, box %g\n", rc, rs, lj->l*.5f);
  }
  if (initload) lj_readpos(lj, lj->x, lj->v, fnpos, LJ_LOADBOX);
  foo(lj);

  md_shiftcom(lj->v, lj->n, lj->d);
  md_shiftang(lj->x, lj->v, lj->n, lj->d);

  for (t = 0; t < nsteps; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    if ((t + 1) % 10 == 0) {
      //vacc += lj_volmove(lj, 1e-2, tp, pressure);
      vacc += lj_lgvvolmove(lj, 3e-4, tp, pressure, 0.2);
      vtot++;
    }
    lj_vrescale(lj, tp, thermdt);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj_calcp(lj, tp));
      if (usesw) {
        bc = lj_bconfsw3d(lj, NULL);
        av_add(&avbc, bc);
      }
    }
    //printf("epot = %g\n", lj->epot); if (lj->epot > 100) exit(1);
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  bc = av_getave(&avbc);
  printf("U/N %6.3f, K/N %6.3f, p %6.3f, bc %6.3f, rho %6.3f, vacc %g%%\n",
      u, k, p, bc, lj->n/lj->vol, 100.0*vacc/vtot);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  return 0;
}
