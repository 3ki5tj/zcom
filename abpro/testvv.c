/* velocity verlet */
#define HAVE_REAL 1
//typedef float real;
typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id = 7, d = 3, model = 2;
  int it, itmax = 1000000;
  real dt = 5e-3f, tp = 0.5f; 
  double smt = 0, sme = 0;

  ab = ab_open(id, d, model);
  ab_initpos(ab, ab->x, 1.0);

  for (it = 1; it <= itmax; it++) {
    ab_vv(ab, 1.f, dt, AB_SOFTFORCE|AB_MILCSHAKE);
    if (it % 10 == 0) ab_rmcom(ab, ab->x, ab->v);
    ab_vrescale(ab, tp, 1e-2);
    if (it % 20000 == 0) {
      printf("t = %4g, T = %6.4f, E = %+8.5f%+8.5f = %+8.5f\n", 
          ab->t, ab->tkin, ab->epot, ab->ekin, ab->epot + ab->ekin);
      //if (it < 30000){ int i; for (i = 0; i < ab->n*ab->d; i++) ab->v[i] *= 0.8f; }
    }
    smt += ab->tkin;
    sme += ab->epot;
  }
  printf("tkin = %g, epot = %g\n", smt/itmax, sme/itmax);
  ab_writepos(ab, ab->x, ab->v, "ab.pos");
  ab_close(ab);
  return 0;
}

