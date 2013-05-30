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

static void copy(real d[3][3], real s[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) d[i][j] = s[i][j];
}

static void rot(real r[3][3], real u[3][3])
{
  real x[3][3];
  rm3_mult(x, r, u);
  copy(r, x);
}

int main(void)
{
  real r[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, u[3][3], v[3][3], detm, amp = 0.5;
  int i;

  for (i = 0; i < 1000; i++) {
    mkrotx(u, amp*rand()/RAND_MAX);
    rot(r, u);
    mkroty(u, amp*rand()/RAND_MAX);
    rot(r, u);
    mkrotz(u, amp*rand()/RAND_MAX);
    rot(r, u);

    detm = rm3_det(r);
    if (fabs(detm - 1) > 0.01) {
      fprintf(stderr, "%d: determinant test failed %g\n", i, detm);
      rm3_print(r, "rot", "%8.3f", 1);
      rm3_print(u, "ro0", "%8.3f", 1);
      break;
    }
    rm3_inv(v, r);
    rm3_mul(u, v, r);
    if (fabs(u[0][0] - 1) > 0.01 || fabs(u[1][1] - 1) > 0.01 || fabs(u[2][2] - 1) > 0.01
     || fabs(u[0][1]) > 0.01 || fabs(u[1][0]) > 0.01
     || fabs(u[1][2]) > 0.01 || fabs(u[2][1]) > 0.01
     || fabs(u[2][0]) > 0.01 || fabs(u[0][2]) > 0.01) {
      fprintf(stderr, "%d: inversion test failed\n", i);
      rm3_print(r, "rot", "%8.3f", 1);
      rm3_print(v, "inv", "%8.3f", 1);
      rm3_print(u, "uni", "%8.3f", 1);
      break;
    }
  }
  printf("passed det/inv tests\n");
  return 0;
}
