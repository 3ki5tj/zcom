#define HAVEREAL 1
typedef float real;
#include "rv3.c"

#define N 128
real max;



/* compute RMSD of two random configurations */
real foo(int n, int nt)
{
  int i, j;
  real x[N][3], xf[N][3], y[N][3], r[3][3], t[3], rmsd1, rmsd2, tmp;

  /* initial configuration */
  for (i = 0; i < n; i++)
    for (j = 0; j < 3; j++) {
      y[i][j] = 1.*rand()/RAND_MAX;
      x[i][j] = 1.*rand()/RAND_MAX;
    }

  rmsd1 = rv3_rmsd(x, NULL, y, NULL, n, r, t);

  for (rmsd2 = 0, i = 0; i < n; i++) {
    real xs[3];
    rv3_add(xf[i], rm3_mulvec(xs, r, x[i]), t); /* xf = R x + t */
    rmsd2 += rv3_dist2(y[i], xf[i]);
  }
  rmsd2 = (real) sqrt(rmsd2/n);
  tmp = fabs(rmsd1 - rmsd2);
  if (tmp > max) max = tmp;
  /* verify rmsd */
  printf("%d, n = %d rmsd %g, corr. %g, diff = %e\n", nt, n, rmsd1, rmsd2, tmp);

  if (fabs(rmsd1 - rmsd2) > 1e-4) {
    if (n == 2) {
      real x01 = rv3_dist(x[0], x[1]), y01 = rv3_dist(y[0], y[1]);
      rv3_print(x[0], "x0", "%22.14f", 1);
      rv3_print(x[1], "x1", "%22.14f", 1);
      rv3_print(y[0], "y0", "%22.14f", 1);
      rv3_print(y[1], "y1", "%22.14f", 1);
      printf("dist01 = %g (x), %g (y), diff %g\n", x01, y01, fabs(x01-y01)/2);
    }
    printf("rmsd mismatch! %g, n %d\n", rmsd1 - rmsd2, n);
    exit(1);
  }
  return rmsd2;
}

int main(void)
{
  int n, t;

  for (t = 1; t <= 1000000; t++) {
    n = 1 + (int)(8.0*rand()/RAND_MAX);
    foo(n, t);
  }
  printf("max rmsd diff = %g\n", max);
  return 0;
}
