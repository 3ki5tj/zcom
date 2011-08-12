/* check constraint algorithms */
#define HAVE_REAL 1
typedef float real;
//typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id, d = 3, model = 1, soft, i, di, j;
  char fn[FILENAME_MAX];
  real Ep, E, E0, E1, del = 0.1f, f2, mag;

  //for (id = 6; id <= 10; id++)
  id = 10;
  {
    ab = ab_open(id, d, model);
    sprintf(fn, "data/%ddm%d/%dbest", d, model, ab->n);
    for (soft = 0; soft < 2; soft++) {
      if (ab_readpos(ab, ab->x, NULL, fn) == 0) {
        Ep = ab_energy(ab, ab->x, soft);
        E = ab_force(ab, ab->f, ab->x, soft);

        /* slightly distub around the equilibrium position */
        for (i = 0; i < ab->n; i++) {
          for (j = 0; j < ab->d; j++) 
            ab->x[i*ab->d + j] += del*(2.f*rand()/RAND_MAX - 1);
        }
        E0 = ab_force(ab, ab->f, ab->x, soft);
        
      }
    }
    ab_close(ab);
  }
  return 0;
}

