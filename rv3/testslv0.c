/* test rv3_solvezero */
#include "rv3.c"

static void foo(real a00, real a01, real a02,
                real a10, real a11, real a12,
                real a20, real a21, real a22)
{
  real a[3][3], as[3][3], x[3][3], y[3][3], res[3];
  int i, n;

  rv3_make(a[0], a00, a01, a02);
  rv3_make(a[1], a10, a11, a12);
  rv3_make(a[2], a20, a21, a22);
  rm3_copy(as, a);
  rm3_print(a, "A", "%10.6f", 1);

  n = rm3_solvezero(a, x, 1e-6);

  printf("%d solutions\n", n);
  for (i = 0; i < n; i++) {
    rm3_mulvec(y[i], as, x[i]);
    rv3_print(x[i], "x", "%10.6f", 1);
    rv3_print(y[i], "y", "%10.6f", 1);
    res[i] = rv3_sqr(y[i]);
    printf("residue %e\n", sqrt(res[i]));
  }
}


int main(void)
{
  foo(  1,  3.,   1.5,
        1,  5.3,  2.5,
        2,  6,    3);
  return 0;
}
