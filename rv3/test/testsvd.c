#include <stdio.h>
#include "rv3.h"

int main(void)
{
  real a[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  real u[3][3], sig[3], v[3][3], r[3][3];
  real detu, detv, detr;

  mat3_svd(a, u, sig, v);
  printf("U:\n%10.5f%10.5f%10.5f\n%10.5f%10.5f%10.5f\n%10.5f%10.5f%10.5f\n",
      u[0][0],u[0][1],u[0][2],u[1][0],u[1][1],u[1][2],u[2][0],u[2][1],u[2][2]);
  printf("V:\n%10.5f%10.5f%10.5f\n%10.5f%10.5f%10.5f\n%10.5f%10.5f%10.5f\n",
      v[0][0],v[0][1],v[0][2],v[1][0],v[1][1],v[1][2],v[2][0],v[2][1],v[2][2]);
  detu = mat3_det(u);
  detv = mat3_det(v);
  mat3_mult(r, v, u);
  detr = mat3_det(r);
  printf("det = %g %g %g\n", detu, detv, detr);
  return 0;
}
