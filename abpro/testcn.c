/* check constraint algorithms */
#define HAVE_REAL 1
typedef float real;
//typedef double real;
#include "abpro.c"

int main(void)
{
  abpro_t *ab;
  int id, d = 3, model = 1, i, j, it, itmax = 100000, milc = 0;
  char fn[FILENAME_MAX];
  real del = 0.1f;

  //for (id = 6; id <= 10; id++)
  id = 9;
  {
    ab = ab_open(id, d, model);
    sprintf(fn, "data/%ddm%d/%dbest", d, model, ab->n);
    die_if (ab_readpos(ab, ab->x, NULL, fn) != 0, "cannot load coordinates");
    memcpy(ab->f, ab->x, ab->n*ab->d*sizeof(real));
    for (it = 0; it < itmax; it++) {
      memcpy(ab->x, ab->f, ab->n*ab->d*sizeof(real));

      /* slightly distub around the equilibrium position */
      for (i = 0; i < ab->n; i++) {
        for (j = 0; j < ab->d; j++) 
          ab->x1[i*ab->d + j] = ab->x[i*ab->d + j] + del*(2.f*rand()/RAND_MAX - 1);
      }
      if (milc) {
        i = ab_milcshake(ab, ab->x, ab->x1, NULL, 0, 0, 0);
      } else {
        i = ab_shake(ab, ab->x, ab->x1, 0, 0);
      }
      die_if (0 != i, "failed to shake properly\n");
      die_if (ab_checkconn(ab, ab->x1, 0) != 0,
        "it %d: broken after a move\n", it);
    }
    printf("completed %d shakes\n", itmax);
    ab_close(ab);
  }
  return 0;
}

