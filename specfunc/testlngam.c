#include <stdlib.h>
#include "specfunc.c"
        
int main(int argc, char **argv)
{
  double x = .5;

  if (argc > 1) x = atof(argv[1]);
  printf("lngamma(%g) = %.14f\n", x, lngam(x));
  return 0;
}

