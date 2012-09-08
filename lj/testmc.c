#include "lj.c"
#include "include/av.h"

int n = 108;
int d = 3;
real rho = 0.8f;
real rcdef = 2.5f;
real tp = 1.0f;
real amp = 0.04f;
int nsteps = 100000;
int usesq = 1;
real ra = 1.0f, rb = 1.25f;
av_t avU, avp;

int main(void)
{
  int t, acc = 0;
  lj_t *lj = lj_open(n, d, rho, rcdef);
  real u, p;

  if (usesq) lj_initsq(lj, ra, rb);
  for (t = 0; t < nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
    if (t >= nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avp, lj_calcp(lj, tp));
    }
  }
  u = av_getave(&avU)/lj->n;
  p = av_getave(&avp);
  printf("erg %g, p %g, acc %g\n", u, p, 1.*acc/nsteps);
  lj_close(lj);
  return 0;
}

