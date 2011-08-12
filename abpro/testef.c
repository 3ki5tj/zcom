/* check energy and force */
#define HAVE_REAL 1
//typedef float real;
typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id, d = 3, model = 2, rndcfg = 0, soft, i, di, j;
  char fn[FILENAME_MAX];
  real Ep, E, E0, E1, del = 0.001f, f2, mag;

  for (id = 6; id <= 10; id++)
  {
    ab = ab_open(id, d, model);
    sprintf(fn, "data/%ddm%d/%dbest", d, model, ab->n);
    for (soft = 0; soft < 2; soft++) {
      if (rndcfg) {
        ab_initpos(ab, ab->x, 1.0);
      } else {
        if (ab_readpos(ab, ab->x, NULL, fn) != 0) {
          fprintf(stderr, "cannot open %s\n", fn);
          continue;
        }
      }
      Ep = ab_energy(ab, ab->x, soft);
      E = ab_force(ab, ab->f, ab->x, soft);

      /* slightly distub around the equilibrium position */
      for (i = 0; i < ab->n; i++) {
        for (j = 0; j < ab->d; j++) 
          ab->x[i*ab->d + j] += del*(2.f*rand()/RAND_MAX - 1);
      }
      E0 = ab_force(ab, ab->f, ab->x, soft);
      
      /* update position according to the force by a unit */
      for (f2 = 0., i = 0; i < ab->n; i++) {
        if (ab->d == 3)
          f2 += rv3_sqr(ab->f + i*ab->d);
        else 
          f2 += rv2_sqr(ab->f + i*ab->d);
      }
      mag = del/f2;
      for (i = 0; i < ab->n; i++) {
        di = i*ab->d;
        if (ab->d == 3)
          rv3_lincomb2(ab->x1 + di, ab->x + di, ab->f + di, 1, mag);
        else
          rv2_lincomb2(ab->x1 + di, ab->x + di, ab->f + di, 1, mag);
      } 
      E1 = ab_energy(ab, ab->x1, soft);
      printf("%dD, model %d, %dmer, soft = %d, E = %.6f, %.6f, f2 = %g, dE = %g\n", 
        d, model, ab->n, soft, Ep, E, f2, (E0 - E1)/del);
    }
    ab_close(ab);
  }
  return 0;
}

