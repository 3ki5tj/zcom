#include <stdio.h>
#include <string.h>

typedef float real;
#define ZCHAVEREAL 1
#include "rv3.h" 

#define R2D (180.0/M_PI)

int main(void)
{
  dihcalc_t dih;

  real a[3] = {1, 0, 0},
       b[3] = {0, 0, 0},
       c[3] = {0, -1, 0},
       d[3] = {-1, -1, 2e-2};

/*
  real a[3] = {1.0, 3.123456789, -0.577350269}, 
       b[3] = {0, 0, 0}, 
       c[3] = {0, -1, 0}, 
       d[3] = {1.0, -1.2345, 1.73205080757};
*/
  real a2[3], b2[3], c2[3], d2[3];
  real g2a, g2b, g2c, g2d, g2;
  real phi, p1, p2, p3, p4, p5;
  real del = 0.001;

  memset(&dih, '\0', sizeof(dih));
  dih.szreal = sizeof (real);
  phi = rv3_calcdih(&dih, a, b, c, d, DIH_GRAD);
  printf("phi = %.14f PI = %.14f\n", phi, (real)(M_PI));
  
  /* test if moving along gradient increases the dihedral by 1.0 */
  g2a = rv3_dot(dih.g[0], dih.g[0]);
  rv3_lincomb2(a2, a, dih.g[0], 1.0, del/g2a);
  p1 = rv3_calcdih(NULL, a2, b, c, d, 0);
  printf("p1 = %g, %g, dphi = %g*del\n", p1, p1*R2D, (p1-phi)/del);

  g2b = rv3_dot(dih.g[1], dih.g[1]);
  rv3_lincomb2(b2, b, dih.g[1], 1.0, del/g2b);
  p2 = rv3_calcdih(NULL, a, b2, c, d, 0);
  printf("p2 = %g, %g, dphi = %g*del\n", p2, p2*R2D, (p2-phi)/del);

  g2c = rv3_dot(dih.g[2], dih.g[2]);
  rv3_lincomb2(c2, c, dih.g[2], 1.0, del/g2c);
  p3 = rv3_calcdih(NULL, a, b, c2, d, 0);
  printf("p3 = %g, %g, dphi = %g*del\n", p3, p3*R2D, (p3-phi)/del);

  g2d = rv3_dot(dih.g[3], dih.g[3]);
  rv3_lincomb2(d2, d, dih.g[3], 1.0, del/g2d);
  p4 = rv3_calcdih(NULL, a, b, c, d2, 0);
  printf("p4 = %g, %g, dphi = %g*del\n", p4, p4*R2D, (p4-phi)/del);

  g2 = g2a + g2b + g2c + g2d;
  rv3_lincomb2(a2, a, dih.g[0], 1.0, del/g2);
  rv3_lincomb2(b2, b, dih.g[1], 1.0, del/g2);
  rv3_lincomb2(c2, c, dih.g[2], 1.0, del/g2);
  rv3_lincomb2(d2, d, dih.g[3], 1.0, del/g2);
  p5 = rv3_calcdih(NULL, a2, b2, c2, d2, 0);
  printf("p5 = %g, %g, dphi = %g*del\n", p5, p5*R2D, (p5-phi)/del);
  return 0;
}

