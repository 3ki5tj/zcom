/* Brownian dynamics */
#define HAVE_REAL 1
typedef float real;
//typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id = 7, d = 3, model = 2;
  int it, itmax = 1000000, n;
  real dt = 1e-3f, T = .5f, sme = 0.;

  ab = ab_open(id, d, model, 0.1);
  n = ab->n;

  for (it = 1; it <= itmax; it++) {
    ab_brownian(ab, T, 1.f, dt, AB_SOFTFORCE|AB_MILCSHAKE);
    if (it % 10000 == 0) 
      printf("t = %4g, e = %+8.5f\n", ab->t, ab->epot);
    sme += ab->epot;
  }
  printf("e = %g\n", sme/itmax);
  ab_writepos(ab, ab->x, NULL, "ab.pos");
  ab_close(ab);
  return 0;
}

