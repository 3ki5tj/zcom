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
  uint32_t x;

  for (i = 0; i < 2*N; i++) {
    x = mtrand();
    id = i % N;
    if (id < 10)
      printf("%6d: %u\n", id, x);
  }
  return 0;
}

