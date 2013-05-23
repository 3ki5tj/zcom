#include <stdio.h>
#include "rng.c"
#define LB 5
#define L  (1<<LB)
#define IX(ix, iy) (((ix)%L) + (iy)*L)
#define IY(ix, iy) ((ix) + ((iy)%L)*L)
#define N (L*L)
int s[N];
double beta = 0.4;

int main(void)
{
  unsigned t;
  int E, h;
  unsigned id = 0, ix, iy;
  uint32_t proba[5] = {0};
  double x, s1 = 0.0, se = 0.0, se2 = 0.0;

  for (id = 0; id < N; id++) s[id] = -1; /* initialize a spin configuration */
  E = -2*N;
  proba[2] = (uint32_t) (4294967295. * exp(-4.*beta));
  proba[4] = (uint32_t) (4294967295. * exp(-8.*beta));
  for (t = 1; t <= 1000000000; t++) {
    /* id = mtrand() >> (32-2*LB); */
    id = (id + 1) % N;
    ix = id % L, iy = id / L; /* randomly choose a spin */
    h = s[id]*(s[IX(ix+1, iy)] + s[IX(ix-1+L, iy)] + s[IY(ix, iy+1)] + s[IY(ix, iy-1+L)]);
    if (h <= 0 || mtrand() <= proba[h]) { /* accept the move */
      E += h * 2;
      s[id] = -s[id];
    }
    if (t % N == 0) {
      s1 += 1.0;
      se += (x = E);
      se2 += x*x;
    }
  }
  x = se / s1;
  printf("beta = %g, E = %d, Eav = %g, Cv = %g\n",
      beta, E, x, (beta*beta)*(se2/s1 - x*x));
  return 0;
}

