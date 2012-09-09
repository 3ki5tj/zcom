#include "lj.c"
#include "include/av.h"

int n = 108;
int d = 3;
real rho = 0.8f;
real rcdef = 2.5f;
real tp = 1.0f;
real amp = 0.04f;
int nsteps = 100000, nstrdf = 100;
int usesq = 1;
real ra = 1.0f, rb = 1.5f;
av_t avU, avp;

int main(void)
{
  int t, acc = 0;
  lj_t *lj = lj_open(n, d, rho, rcdef);
  real u, p;

  if (usesq) lj_initsq(lj, ra, rb);
  lj_rdfopen(lj, 0.01);
  for (t = 1; t <= nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avp, lj_calcp(lj, tp));
      if (t % nstrdf == 0) lj_rdfadd(lj);
    }
  }
  u = av_getave(&avU)/lj->n;
  p = av_getave(&avp);
  printf("erg %g, p %g, acc %g, rdfnfr %d\n", u, p, 1.*acc/nsteps, lj->rdfnfr);
  lj_rdfsave(lj, "rdf.dat", HIST_NOZEROES);
  lj_close(lj);
  return 0;
}

