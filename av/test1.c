#include "av.h"
#define N 8
int main(void)
{
  double arr[N] = {0, -0.4, .3, 0, 0, 0, -.3, .4};
  int i;
  av_t s[1];

  av_clear(s);
  for (i = 0; i < N; i++) av_add(s, arr[i]);
  printf("average %g, standard dev. %g\n", av_getave(s), av_getdev(s));
  return 0;
}
