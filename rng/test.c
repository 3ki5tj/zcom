#include <stdio.h>

#ifdef REF
#include "zcomref.h"
#else
#include "rng.c"
#endif

#ifndef N
#define N 1000000
#endif

int main(void) 
{
  int i, id;
  UINT32 x;

  for (i = 0; i < 2*N; i++) {
    x = mtrand();
    id = i % N;
    if (id < 10)
      printf("%6d: " UI32FMT "\n", id, x);
  }
  return 0;
}

