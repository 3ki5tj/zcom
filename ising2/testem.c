/* check if pick() and flip() correctly computes the energy */
#define L 32
/*
#define IS2_LB 5
*/
#include "ising2.h"

/* randomly pick a site and flip it */
static void randflip(ising_t *is)
{
  int id, h, E, M, E2, M2;
/*
  id = is2_pick(is, &h);
  is2_flip(is, id, h);
*/
  IS2_PICK(is, id, h);
  IS2_FLIP(is, id, h);
  E = is->E;
  M = is->M;
  E2 = is2_em(is);
  M2 = is->M;
  printf("id: %4d, s = %d, h = %2d, E = %4d, %4d, M = %4d, %4d\n",
      id, is->s[id], h, E, E2, M, M2);
  if (is2_check(is) != 0) exit(1);
}

int main(void)
{
  ising_t *is;
  int i;

  if ((is = is2_open(L)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  for (i = 0; i < 10; i++) randflip(is);
  printf("E = %d, M = %d\n", is->E, is->M);
  is2_close(is);
  return 0;
}

