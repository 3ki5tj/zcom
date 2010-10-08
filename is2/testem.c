#include <stdio.h>
#include "is2.h"

/* randomly pick a site and flip it */
static void randflip(is_t *is)
{
  int id, nd, E, M, E2, M2;

  id = is2_pick(is, &nd);
  E = is2_flip(is, id, nd);
  M = is->M;
  E2 = is2_em(is);
  M2 = is->M;
  printf("id: %4d, nd = %2d, E = %4d, %4d, M = %4d, %4d\n", 
      id, nd, E, E2, M, M2);
  if (E != E2 || M != M2) {
    fprintf(stderr, "error occured\n");
    exit(1);
  }
}

int main(void)
{
  is_t *is;
  int i;

  if ((is = is2_open(5)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }

  for (i = 0; i < 10; i++) randflip(is);
  
  printf("E = %d, M = %d\n", is->E, is->M);
  is2_close(is); 
  return 0;
}

