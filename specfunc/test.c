#include <stdlib.h>
#include "specfunc.c"

int main(int argc, char **argv)
{
  double x = .5;
  if (argc > 1) x = atof(argv[1]);
  printf("lngamma(%g) = %.14f\n", x, lngam(x));
  printf("sqrt(2pi) = %.18f\n", sqrt(2*M_PI));
  return 0;
}

