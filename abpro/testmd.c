/* molecular dynamics */
#define HAVE_REAL 1
/*
typedef float real;
*/
typedef double real;
#include "abpro.c"

int doconstr = 1;
int br = 0;

int main(void)
{
  abpro_t *ab;
  int id = 10, d = 3, model = 1;
  int it, itmax = 100000;
  real brdt = 1e-4f, mddt = 2e-3f, tp = 0.4f, thermdt = 0.01f; 
  double smt = 0, sme = 0;
  unsigned flags = AB_SOFTFORCE;
  
  ab = ab_open(id, d, model, 0.1);
  if (doconstr) { /* optimal constraints */
    ab_initconstr(ab, -1);
  } else {
    flags |= AB_MILCSHAKE;
  }
 
  for (it = 1; it <= itmax; it++) {
    if (br) {
      ab_brownian(ab, tp, 1.f, brdt, flags);
    } else {
      ab_vv(ab, 1.f, mddt, flags);
      ab_rmcom(ab, ab->x, ab->v);
      ab_vrescale(ab, tp, thermdt);
    }
    if (it % 1000 == 0 && ab->lgcon && ab->lgact < ab->lgcnt)
      ab_updconstr(ab, 0);
    if (it % 10000 == 0) 
    {
      printf("t = %4g, T = %6.4f, E = %+9.4f%+9.4f = %+9.4f constr %d/%d\n", 
          ab->t, ab->tkin, ab->epot, ab->ekin, ab->epot + ab->ekin, ab->lgact, ab->lgcnt);
    }
    smt += ab->tkin;
    sme += ab->epot;
  }
  printf("tkin = %g, epot = %g\n", smt/itmax, sme/itmax);
  //ab_writepos(ab, ab->x, br ? NULL : ab->v, "ab.pos");
  //ab_close(ab);
  return 0;
}

