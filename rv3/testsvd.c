/*
#define HAVEREAL 1
typedef float real;
*/
#define RV3_DEBUG
#include "rv3.c"
#include "svd.h"

/*
//real a[3][3] = {{1, 2, 0},{2, 4, 0}, {0, 0, 5}};
//real a[3][3] = {{1, 2, 0},{2, 4, 0}, {0, 0, 5+1e-8}};
//real a[3][3] = {{1, 2, 0},{2, 4, 0}, {0, 0, 5+1e-9}};
//real a[3][3] = {{1, 2, 0},{2, 4, 0}, {0, 0, 5-1e-8}};
*/


real a[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};


/*
real a[3][3] = {
{   -0.05425476794102,    0.05788731990739,    0.00650199809425},
{    0.15857706314977,   -0.16919436821669,   -0.01900418712532},
{    0.06902557678025,   -0.07364708755577,   -0.00827216087566 }};
*/

real u[3][3], sig[3], v[3][3], r[3][3];
real detu, detv, detr;

static void dump(void)
{
  const char *rfmt = "%16.10f";
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
  rv3_print(sig, "S", "%20.14f", 1);
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
