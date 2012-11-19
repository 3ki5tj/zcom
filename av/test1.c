#include "av.h"
#define N 8
int main(void)
{
  double arr[N] = {0, -0.4, .3, 0, 0, 0, -.3, .4};
  double arr2[N] = {100, 110, 125, 90, 75, 105, 85, 110};
  int i;
  av_t s[1];
  avn_t *avn;

  av_clear(s);
  for (i = 0; i < N; i++) av_add(s, arr[i]);
  printf("average %g, standard dev. %g\n", av_getave(s), av_getdev(s));

  avn = avn_open(2);
  for (i = 0; i < N; i++) {
    avn_addw(avn, 1, arr[i], arr2[i]);
  }
  printf("arr1: ave %g, dev. %g\n", avn_getave(avn, 0), avn_getdev(avn, 0));
  printf("arr2: ave %g, dev. %g\n", avn_getave(avn, 1), avn_getdev(avn, 1));

  return 0;
}
