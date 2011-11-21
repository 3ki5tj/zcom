/* check if pick() and flip() correctly computes the energy */
#include <stdio.h>
#include "rng.c"


#define PT2_LB 5
//#define PT2_L 32
#define PT2_Q 10
#include "potts2.c"

/* randomly pick a site and flip it */
static void randflip(potts_t *pt)
{
  int id, so, sn, h[PT2_Q], E, M, E2, M2;
/*
  id = pt2_pick(pt, h);
  PT2_NEWFACE(pt, id, so, sn);
  pt2_flip(pt, id, sn, h);
*/

  PT2_PICK(pt, id, h);
  PT2_NEWFACE(pt, id, so, sn);
  PT2_FLIP(pt, id, so, sn, h);

  E = pt->E;
  M = pt->M[0];
  E2 = pt2_em(pt);
  M2 = pt->M[0];
  printf("id: %4d, s = %d, h = %2d, E = %4d, %4d, M = %4d, %4d\n", 
      id, pt->s[id], h[0], E, E2, M, M2);
  if (pt2_check(pt) != 0) exit(1);
}

int main(void)
{
  potts_t *pt;
  int i;

  if ((pt = pt2_open(PT2_L, PT2_Q)) == NULL) {
    fprintf(stderr, "cannot init\n");
    return -1;
  }
  for (i = 0; i < 1000; i++) randflip(pt);
  printf("E = %d, M = %d\n", pt->E, pt->M[0]);
  pt2_close(pt); 
  return 0;
}

