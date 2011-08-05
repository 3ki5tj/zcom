#include <stdio.h>
#include "rng.c"
#include "ising2.h"

int main(void)
{
  ising_t *is;
  double beta = 1e-6, eav, cv, lnz;

  if ((is = is2_open(5)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  lnz = is2_exact(is, beta, &eav, &cv);
  printf("lnZ: %.15f, Eav: %.15e, Cv: %.15e\n", lnz, eav, cv);
  is2_close(is); 
  return 0;
}

