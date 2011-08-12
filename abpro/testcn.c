/* check constraint algorithms */
#define HAVE_REAL 1
//typedef float real;
typedef double real;
#include "abpro.c"


/* randomly displace ab->x by del, then shake it
 * save result to ab->x1 (if repl, also to ab->x) */
static int randmove(abpro_t *ab, real del, int repl, 
    int milc, int it, const char *tag)
{
  int i, j, verbose = 0;

  /* slightly distub around the equilibrium position */
  for (i = 0; i < ab->n; i++) {
    for (j = 0; j < ab->d; j++) 
      ab->x1[i*ab->d + j] = ab->x[i*ab->d + j] + del*(2.f*rand()/RAND_MAX - 1);
  }
  if (milc) {
    i = ab_milcshake(ab, ab->x, ab->x1, NULL, 0.f, 0, 0, verbose);
  } else {
    i = ab_shake(ab, ab->x, ab->x1, 0, 0.f, verbose);
  }
  die_if (0 != i, "%s| it %d: failed to shake properly\n", tag, it);
  die_if (ab_checkconn(ab, ab->x1, 0) != 0,
    "%s| it %d: broken after a move\n", tag, it);
  if (repl) {
    memcpy(ab->x, ab->x1, ab->n*ab->d*sizeof(real));
  }
  return 0;
}

int main(void)
{
  abpro_t *ab;
  int id, d = 2, model = 1, it, itmax = 100000, milc = 1, isrand = 1;
  char fn[FILENAME_MAX];
  real del = 0.1f;

  //for (id = 6; id <= 10; id++)
  id = 10;
  {
    ab = ab_open(id, d, model);
    if (isrand) {
      ab_initpos(ab, ab->x, 5.0);
    } else { /* save init structure to ab->xmin */
      sprintf(fn, "data/%ddm%d/%dbest", d, model, ab->n);
      if (ab_readpos(ab, ab->x, NULL, fn) != 0) {
        fprintf(stderr, "Warning cannot open %s\n", fn);
      }
      memcpy(ab->xmin, ab->x, ab->n*ab->d*sizeof(real));
    }
    for (it = 0; it < itmax; it++) {
      randmove(ab, del, isrand, milc, it, "run");
      if (!isrand) {
        memcpy(ab->x, ab->xmin, ab->n*ab->d*sizeof(real));
      }
    }
    ab_writepos(ab, ab->x, NULL, "ab.pos");
    printf("completed %d shakes\n", itmax);
    ab_close(ab);
  }
  return 0;
}

