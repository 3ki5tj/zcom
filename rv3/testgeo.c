/* test geometric routines
 * rv3_ang
 * rv3_vdist
 * rv3_vpdist
 * */
#include <stdio.h>
#include "rv3.h"


static void testang(void)
{
  real a[3] = {1, 1, 0}, b[3] = {0, 0, 0}, c[3] = {0, 1, 1};
  real fa[3], fb[3], fc[3];
  real ang, ang2, del = .001f;

  ang = rv3_ang(a, b, c, fa, fb, fc);
  printf("ang = %g\n", ang);

  rv3_sinc(a, fa, del/rv3_sqr(fa));
  printf("%g, %g, %g\n", a[0], a[1], a[2]);
  ang2 = rv3_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang2, ang2 - ang);

  rv3_sinc(b, fb, del/rv3_sqr(fb));
  ang = rv3_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang, ang - ang2);

  rv3_sinc(c, fc, del/rv3_sqr(fc));
  ang2 = rv3_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang, ang2 - ang);
}

static void testvdist(void)
{
  real a[3] = {0, 0, 0}, b[3] = {1, 0, 0}, c[3] = {2, 1, 1.732};
  printf("vertical distance from c to a and b is %g\n", rv3_vdist(c, a, b));
}

static void testvpdist(void)
{
  real a[3] = {0, 0, 0},
       b[3] = {1, 0, 0},
       c[3] = {1, 1, 0},
       d[3] = {3, 4, 0.5};
  real vd1, vd2, phi;

  /* compute the distance from d to plane a-b-c */
  /* method 1 */
  vd1 = rv3_vpdist(d, a, b, c);

  /* method 2 */
  vd2 = rv3_vdist(d, a, b);
  phi = rv3_calcdih(NULL, d, a, b, c, 0);
  vd2 *= -sin(phi);

  printf("vd %g, %g; dihedral d - (a-b) - c = %g\n", vd1, vd2, phi);
}

int main(void)
{
  testang();
  testvdist();
  testvpdist();
  return 0;
}
