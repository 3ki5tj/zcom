/* velocity verlet */
#define HAVE_REAL 1
//typedef float real;
typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id = 10, d = 3, model = 2;
  int it, itmax = 100000;
  real dt = 2e-3f, tp = 0.5f, thermdt = 0.01f; 
  double smt = 0, sme = 0;

  ab = ab_open(id, d, model, 0.1);

  for (it = 1; it <= itmax; it++) {
    ab_vv(ab, 1.f, dt, AB_SOFTFORCE|AB_MILCSHAKE);
    if (it % 10 == 0) ab_rmcom(ab, ab->x, ab->v);
    ab_vrescale(ab, tp, thermdt);
    if (it % 20000 == 0) {
      printf("t = %4g, T = %6.4f, E = %+8.5f%+8.5f = %+8.5f\n", 
          ab->t, ab->tkin, ab->epot, ab->ekin, ab->epot + ab->ekin);
    }
    smt += ab->tkin;
    sme += ab->epot;
  }
  printf("tkin = %g, epot = %g\n", smt/itmax, sme/itmax);
  ab_writepos(ab, ab->x, ab->v, "ab.pos");
  ab_close(ab);
  return 0;
}

