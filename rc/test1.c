#include <stdio.h>
#include <math.h>
#include "rc.h"

int main(void)
{
  rcomplex_t a, b, c, d;
  a = rc_make(2, 1);
  b = rc_make(3, 4);
  c = rc_mul(a, b);
  d = rc_div(b, a);

  printf("c = %g + %g i\n", c.re, c.im);
  printf("d = %g + %g i\n", d.re, d.im);
  return 0;
}
