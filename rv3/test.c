#include <stdio.h>
#include "rv3.h"

int main(void)
{
  real a[3] = {0, 0, 0}, b[3] = {1, 0, 0}, c[3] = {2, 1, 1.732};
  printf("vertical distance from c to a and b is %g\n", rv3_vdist(c, a, b));
  return 0;
}
