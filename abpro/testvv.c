/* velocity verlet */
#define HAVE_REAL 1
typedef float real;
//typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id = 6, d = 3, model = 2, milc = 0, soft = 1;
  int it, itmax = 1000000, n;
  real dt = 5e-3f; 

  ab = ab_open(id, d, model);
  ab_initpos(ab, ab->x, 1.0);
  n = ab->n;

  for (it = 1; it <= itmax; it++) {
    ab_vv(ab, dt, soft, milc);
    if (it % 10000 == 0) 
      printf("t = %4g, e = %+8.5f%+8.5f = %+8.5f\n", ab->t, ab->epot/n, ab->ekin/n, (ab->epot + ab->ekin)/n);
  }
  ab_writepos(ab, ab->x, ab->v, "ab.pos");
  ab_close(ab);
  return 0;
}

