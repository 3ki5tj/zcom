#include <stdlib.h>
#include "specfunc.c"
        
static void test_ksq(double x)
{
  FILE *fp;
  printf("ksq(%g) = %.14f\n", x, ksq(x));
  printf("%.18f\n", sqrt(1.23370055013616983));

  /* output KS-distribution for plot, several values:
   * q(0.83) = .5, q(1.22) = .1, q(1.62) = .01 */
  if ((fp = fopen("ksdist", "w")) != NULL) {
    for (x = 0.0; x < 2.0; x += 0.001) {
      double y = ksq(x);
      fprintf(fp, "%g %g %g\n", x, 1-y, y);
    }
    fclose(fp);
  }
} 

int main(int argc, char **argv)
{
  double x = .5;

  if (argc > 1) x = atof(argv[1]);
  test_ksq(x);
  return 0;
}

