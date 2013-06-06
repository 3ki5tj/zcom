#include <stdio.h>
#include "rng.h"

#define N 1000000

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

