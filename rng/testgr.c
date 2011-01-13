#include "rng.c"

#define TMAX 100000000

#define XMIN (-5.0)
#define XMAX (5.0)
#define XDEL  0.01
#define XCNT (int)((XMAX-XMIN)/XDEL + 1)
double hist[XCNT];

int main(void)
{
  int t, i;
  double x;

  for (t = 0; t < TMAX; t++) {
    x = grand0();
    i = (int)((x - XMIN)/XDEL);
    if (i >= 0 && i < XCNT)
      hist[i] += 1.0;
  }
  for (i = 0; i < XCNT; i++)
    printf("%g %g\n", XMIN+XDEL*i, hist[i]/(TMAX*XDEL));
  return 0;
}
