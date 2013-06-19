
#define HAVEREAL 1
typedef float real;

#include "rv3.c"

#define N 5
int main(void)
{
  int i;
  real x[N][3], y[N][3], tmp[3], dev;
  real r0[3][3] = {{1, 0, 0}, {0, 0.8, 0.6}, {0, -.6, .8}},
       t0[3] = {0.75, .25, .5};
  real r[3][3], t[3];

  /* initial configuration */
  for (i = 0; i < N; i++) rv3_rnd0(y[i]);

  /* construct a configuration */
  for (i = 0; i < N; i++) {
    rv3_diff(tmp, y[i], t0);
    rm3_tmulvec(x[i], r0, tmp);
  }

#define DUMP() \
  rm3_print(r,  "R ", "%10.5f", 1); \
  rm3_print(r0, "R0", "%10.5f", 1); \
  rv3_print(t,  "t ", "%10.5f", 1); \
  rv3_print(t0, "t0", "%10.5f", 1);

  printf("\n\n1. Fitting a rotated object using only rotation and translation\n");
  dev = rv3_rmsd(x, NULL, y, NULL, N, 0, r, t);
  printf("dev = %g\n", dev);
  DUMP()

  /* construct a configuration */
  for (i = 0; i < N; i++) {
    rv3_diff(tmp, y[i], t0);
    rm3_tmulvec(x[i], r0, tmp);
    x[i][2] = -x[i][2];
  }

  printf("\n\n2. Fitting a rotated and reflected object using only rotation and translation\n");
  dev = rv3_rmsd(x, NULL, y, NULL, N, 0, r, t);
  printf("dev = %g\n", dev);
  DUMP()

  printf("\n\n3. Fitting a rotated and reflected object using rotation, translation and reflection\n");
  dev = rv3_rmsd(x, NULL, y, NULL, N, 1, r, t);
  printf("dev = %g\n", dev);
  DUMP()

  return 0;
}
