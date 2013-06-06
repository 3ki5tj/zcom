#include <stdio.h>

#include "rng.h"

#define N 10000

int main(void)
{
  int i;
  double v[3];

  for (i = 0; i < N; i++) {
    randor3d(v);
    printf("%g, %g, %g\n", v[0], v[1], v[2]);
  }
  return 0;
}

