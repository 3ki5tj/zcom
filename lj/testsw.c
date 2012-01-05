/* test switched potential */
#define HAVE_REAL 1
typedef float real;
#include "lj.c"

int n = 256;
int d = 3;
real rho = 0.8f;
real rs = 2.0f, rc = 3.0f;
real tp = 1.0f;
real amp = 0.04f;
int nsteps = 100000;

static void printswitch(lj_t *lj, real dr)
{
  real r, u, fscal, psi, xi, lap;

  for (r = dr; r <= lj->rc; r += dr) {
    u = lj_potsw(lj, r, &fscal, &psi, &xi);
    lap = psi * r * r - 3 * fscal;
    printf("%8.4f %.14e %.14e %.14e %.14e %.14e\n",
        r, u, fscal, psi, xi, lap);
  }
}

int main(void)
{
  lj_t *lj = lj_open(n, d, rho, rc);
  lj_initsw(lj, rs);
  printswitch(lj, 0.01f);
  lj_close(lj);
  return 0;
}
