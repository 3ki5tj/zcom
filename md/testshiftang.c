#include "md.h"

#define N 10

static real angmom2d(rv2_t *x, rv2_t *v, int n)
{
  real p = 0, xc[2] = {0, 0}, xi[2];
  int i;

  for ( i = 0; i < n; i++ ) {
    rv2_inc(xc, x[i]);
  }
  rv2_smul(xc, 1.0/n);
  for ( i = 0; i < n; i++ ) {
    rv2_diff(xi, x[i], xc);
    p += xi[0] * v[i][1] - xi[1] * v[i][0];
  }
  return p;
}


int main(void)
{
  rv2_t x[N] = {{0, 0}}, v[N] = {{0, 0}};
  int i;

  for ( i = 0; i < N; i++ ) {
    rv2_rnd0(x[i]);
    rv2_rnd0(v[i]);
  }

  md_shiftcom2d(v, N);
  md_shiftang2d(x, v, N);
  printf("angular momentum, %g\n", angmom2d(x, v, N));
  return 0;
}
