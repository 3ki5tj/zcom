#include <stdio.h>
#include <time.h>
#include "rng.h"



static void testspeed(int nsteps)
{
  clock_t t0 = clock();
  double tsum;
  int i;
  unsigned q = 0;

  mtrng_t *mr = mtrng_open(3254);
  for (i = 0; i < nsteps; i++) {
    q += mtrng_rand(mr);
  }
  tsum = 1. * (clock() - t0) / CLOCKS_PER_SEC;
  printf("mtrand speed: %gs/%d = %gns\n", tsum, nsteps, 1e9*tsum/nsteps);
}



int main(int argc, char **argv)
{
  int nsteps = 2000000000;

  if (argc > 1) nsteps = atoi(argv[1]);
  testspeed(nsteps);
  return 0;
}

