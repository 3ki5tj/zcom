#include <stdio.h>
#include "rv3.h"

int main(void)
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
  
  return 0;
}
