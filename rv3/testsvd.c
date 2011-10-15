#include <stdio.h>
#define HAVEREAL 1
typedef float real;
#include "rv3.h"
#include "svd.c"
  

real a[3][3] = {{1, 2, 3},{4, 5, 6}, {7, 8, 9}};

/*
real a[3][3] = {{1, -1, 0}, {-1, 1, 0}, {0, 0, 2}};

real a[3][3] = {{1.463, -5.472, -14.393}, {-29.442, 27.743, 4.775}, {-77.255, 6.995, 17.999}};
*/
real u[3][3], sig[3], v[3][3], r[3][3];
real detu, detv, detr;

static void dump(void)
{
  const char *rfmt = "%12.6f";
  real usvt[3][3];
  int i;

  rm3_print(a, "A", rfmt, 1);
  rm3_print(u, "U", rfmt, 1);
  rm3_print(v, "V", rfmt, 1);
  detu = rm3_det(u);
  detv = rm3_det(v);
  rm3_mult(r, v, u);
  rm3_print(r, "R", rfmt, 1);
  detr = rm3_det(r);
  rm3_trans(v);
  for (i = 0; i < 3; i++) rv3_smul(v[i], sig[i]);
  rm3_mul(usvt, u, v);
  rm3_print(usvt, "USV^T", rfmt, 1);
  printf("det: u %g, v %g, r %g\n", detu, detv, detr);
  rv3_print(sig, "S", rfmt, 1);
}

int main(void)
{
  printf("Test rm3_svd\n");
  rm3_svd(a, u, sig, v);
  dump();

  printf("\n\nTest svd\n");
  memcpy(u, a, 9*sizeof(real));
  svd((real *) u, sig, (real *) v, 3, 3);
  dump();
  return 0;
}
