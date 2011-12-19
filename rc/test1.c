#include <stdio.h>
#include <math.h>
#include "rc.h"

int main(void)
{
  rcomplex_t a, b, c;
  a = rc_make(1, 2);
  b = rc_make(3, 4);
  c = rc_mul(a, b);
  printf("c = %g + %g i\n", c.re, c.im);
  return 0;
}
