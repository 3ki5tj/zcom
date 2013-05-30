/* test random matrix */
#include <stdio.h>
#include <stdlib.h>
#include "rv3.h"

static void mkrotx(real m[3][3], real th)
{
  int i, j;
  real c = (real) cos(th), s = (real) sin(th);
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) m[i][j] = (i == j);
  m[1][1] =  c; m[1][2] = -s;
  m[2][1] =  s; m[2][2] =  c;
}

static void mkroty(real m[3][3], real th)
{
  int i, j;
  real c = (real) cos(th), s = (real) sin(th);
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) m[i][j] = (i == j);
  m[2][2] =  c; m[2][0] = -s;
  m[0][2] =  s; m[0][0] =  c;
}

static void mkrotz(real m[3][3], real th)
{
  int i, j;
  real c = (real) cos(th), s = (real) sin(th);
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) m[i][j] = (i == j);
  m[0][0] =  c; m[0][1] = -s;
  m[1][0] =  s; m[1][1] =  c;
}

int main(void)
{
  real u[3][3], v[3][3], amp = 0.5f, tol = 1e-5f;
  real vx[3] = {2, 0, 0}, vy[3] = {0, -2.3, 0}, vz[3] = {0, 0, 1.5f}, vd[3] = {1, 1, 1}, a[3];

  mkrotx(u, amp);
  rm3_mkrot(v, vx, amp);
#define check(u, v) \
  if (fabs(u[0][0] - v[0][0]) > tol || fabs(u[0][1] - v[0][1]) > tol || fabs(u[0][2] - v[0][2]) > tol     \
   || fabs(u[1][0] - v[1][0]) > tol || fabs(u[1][1] - v[1][1]) > tol || fabs(u[1][2] - v[1][2]) > tol     \
   || fabs(u[2][0] - v[2][0]) > tol || fabs(u[2][1] - v[2][1]) > tol || fabs(u[2][2] - v[2][2]) > tol) {  \
    printf("u and v mismatch\n"); \
    rm3_print(u, "u", "%8.3f", 1); \
    rm3_print(v, "v", "%8.3f", 1); \
    exit(1); }
  check(u, v);

  mkroty(u, -amp);
  rm3_mkrot(v, vy, amp);
  check(u, v);

  mkrotz(u, amp);
  rm3_mkrot(v, vz, amp);
  check(u, v);

  rv3_print(vx, "x --> y", "%8.3f", 1);
  rv3_rot(a, vx, vd, M_PI*2/3);
  rv3_print(a, "x --> y", "%8.3f", 1);
  return 0;
}
