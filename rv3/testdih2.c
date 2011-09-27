#include "rv3.h"

int main(void)
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
  
  printf("vd %g, %g; phi = %g\n", vd1, vd2, phi);
  return 0;
}

