#include <stdio.h>
#include "rv2.h"

int main(void)
{
  real a[2] = {1, 1.732}, b[2] = {0, 0}, c[2] = {0, 1};
  real fa[2], fb[2], fc[2];
  real ang, ang2, del = .001f;

  ang = rv2_ang(a, b, c, fa, fb, fc);
  printf("ang = %g\n", ang);

  rv2_sinc(a, fa, del/rv2_sqr(fa));
  ang2 = rv2_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang2, ang2 - ang); 
  
  rv2_sinc(b, fb, del/rv2_sqr(fb));
  ang = rv2_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang, ang - ang2); 

  rv2_sinc(c, fc, del/rv2_sqr(fc));
  ang2 = rv2_ang(a, b, c, NULL, NULL, NULL);
  printf("ang = %g, del = %g\n", ang, ang2 - ang); 
  
  return 0;
}
