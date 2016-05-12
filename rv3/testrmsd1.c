
#define HAVEREAL 1
typedef float real;

#include "rv3.c"

#define N 5
int main(void)
{
  int i;
  real x[N][3], y[N][3], yc[3] = {0, 0, 0}, dev;
  real r0[3][3] = {{1, 0, 0}, {0, 0.8, 0.6}, {0, -.6, .8}},
       t0[3] = {0.75, .25, .5};
  real r[3][3], t[3];

  /* initial configuration, y */
  for (i = 0; i < N; i++) {
    rv3_rnd0(y[i]);
    rv3_inc(yc, y[i]);
  }
  rv3_smul(yc, (real) 1/N);
  for ( i = 0; i < N; i++ ) {
    rv3_dec(y[i], yc);
  }

  /* make x by rotating y and translation */
  for (i = 0; i < N; i++) {
    rm3_tmulvec(x[i], r0, y[i]);
    rv3_dec(x[i], t0);
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

  /* make x by rotating y and translation and reflection */
  for (i = 0; i < N; i++) {
    rm3_tmulvec(x[i], r0, y[i]);
    x[i][2] = -x[i][2];
    rv3_dec(x[i], t0);
  }

  printf("\n\n2. Fitting a rotated and _reflected_ object using only rotation and translation\n");
  dev = rv3_rmsd(x, NULL, y, NULL, N, 0, r, t);
  printf("dev = %g\n", dev);
  DUMP()

  printf("\n\n3. Fitting a rotated and _reflected_ object using rotation, translation and reflection\n");
  dev = rv3_rmsd(x, NULL, y, NULL, N, 1, r, t);
  printf("dev = %g\n", dev);
  DUMP()

  return 0;
}
