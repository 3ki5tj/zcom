/* check if the force matches energy */
#include "lj.h"



int N = 108;
real rho = 0.7f;
real rcdef = 1000.0f; /* use half-box cut-off */
int usesw = 0;
real rs = 2.0f;
const char *fnpos = "lj.pos";



/* see if the force matches energy */
static void checkforce(lj_t *lj)
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
  lj_t *lj = lj_open(N, 3, rho, rcdef);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g, box %g\n", lj->rc, rs, lj->l*.5f);
  }
  if (0 == lj_readpos(lj, lj->x, lj->v, fnpos, LJ_LOADBOX))
    fprintf(stderr, "loaded previous coordinates from %s\n", fnpos);
  checkforce(lj);
  lj_close(lj);
  return 0;
}
