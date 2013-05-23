#include <stdlib.h>
#include "specfunc.c"

int main(int argc, char **argv)
{
  double a = .5, x, y;

  if (argc > 1) a = atof(argv[1]);
  printf("lngamma(%g) = %.14f\n", a, lngam(a));
  for (x = 0; x < a; x += .1) {
    y = lnincgam0(a, x);
    printf("lnincgam0(%g, %g) = %g\n", a, x, y);
  }
  for (x = a; x < a + 10; x += 1.) {
    y = lnincgam1(a, x);
    printf("lnincgam1(%g, %g) = %g\n", a, x, y);
  }
  printf("%g\n", lnincgam(159, 200));
  return 0;
}

