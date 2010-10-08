#include <stdio.h>

#ifdef REF
#include "zcomref.h"
#else
#include "rng.h"
#endif

#ifndef N
#define N 1000000
#endif

int main(void) 
{
  int i, id;
  double x;

  for (i = 0; i < 2*N; i++) {
    x = rnd0();
    id = i % N;
    if (id < 10)
      printf("%6d: %16.14f\n", id, x);
  }
  return 0;
}
