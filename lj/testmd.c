#include "lj.c"
#include "include/av.h"

int N = 256;
real rho = 0.8f;
real rc = 3.0f;
real tp = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;
int usesw = 0;
real rs = 2.0f;
av_t avU, avK, avp;

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
  int t, nsteps = 20000;
  real u, k, p;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g\n", rc, rs);
  }
  foo(lj);
  for (t = 0; t < nsteps; t++) {
    lj_vv(lj, mddt);
    lj_vrescale(lj, tp, thermdt);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj->rho * tp + lj->pvir);
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp)/N;
  printf("U/N = %6.3f, K/N = %6.3f, p = %6.3f\n", u, k, p);
  lj_close(lj);
  return 0;
}
