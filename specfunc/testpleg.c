#include "specfunc.h"
int main(void)
{
  double x, c, y20, y22, y1, y2;

  for (x = 0; x < M_PI; x += 0.1) {
    c = cos(x);
    y20 = sqrt(5./4/M_PI)*.5*(c*c*3-1);
    y22 = sqrt(45./4/M_PI/24)*(1-c*c);
    y1 = plegendre(c, 2, 0);
    y2 = plegendre(c, 2, 2);
    printf("%g %g %g %g %g\n", x, y1, y20, y2, y22);
  }
  return 0;
}
