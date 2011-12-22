#include <stdio.h>
#include <math.h>
#include "distr.c"

#define N 360
#define M 360
#define N1 (N+1)
#define M1 (M+1)
real dx = 2*M_PI/N, dy = 2*M_PI/M;

int main(void)
{
  distr2d_t *d = distr2d_open(-M_PI, M_PI, dx, -M_PI, M_PI, dy);
  double x, y;
  int i, j, id1;

  for (i = 0; i < N; i++) {
    x = d->xmin + (i + .5) * dx;
    for (j = 0; j < M; j++) {
      y = d->ymin + (j + .5) * dy;
      id1 = i*M1 + j;
      d->mf[id1] = x;
      d->mg[id1] = 2*y;
    }
  }
 
  distr2d_calclnrho(d);
  distr2d_save(d, "ab.txt");
  distr2d_close(d); 
  return 0;
}
