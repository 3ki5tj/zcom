/* Brownian dynamics */
#define HAVE_REAL 1
typedef float real;
//typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id = 6, d = 3, model = 2, milc = 0, soft = 1;
  int it, itmax = 1000000, n;
  real dt = 5e-4f, T = 0.3; 

  ab = ab_open(id, d, model);
  ab_initpos(ab, ab->x, 1.0);
  n = ab->n;

  for (it = 1; it <= itmax; it++) {
    ab_brownian(ab, T, dt, soft, milc);
    if (it % 10000 == 0) 
      printf("t = %4g, e = %+8.5f\n", ab->t, ab->epot/n);
  }
  ab_writepos(ab, ab->x, NULL, "ab.pos");
  ab_close(ab);
  return 0;
}

