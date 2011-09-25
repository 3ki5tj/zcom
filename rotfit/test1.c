#include <stdlib.h>
#include "rotfit.c"

#define N 5
int main(void)
{
  int i, j;
  real x[N][3], y[N][3], tmp[3], dev;
  real r0[3][3] = {{1, 0, 0}, {0, 0.8, 0.6}, {0, -.6, .8}}, 
       t0[3] = {0.75, .25, .5};
  real r[3][3], t[3];

  /* initial configuration */
  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {
      y[i][j] = 1.*rand()/RAND_MAX;
    }
  }

  /* construct a configuration */
  for (i = 0; i < N; i++) {
    rv3_diff(tmp, y[i], t0);
    rm3_mulvec(x[i], r0, tmp);
  }
  
  printf("\n\nPure rotation\n");
  dev = rotfit3(x, NULL, y, NULL, N, r, t);  
  printf("dev = %g\n", dev);
  rm3_print(r, "R ", "%10.5f", 1);
  rm3_print(r0, "R0", "%10.5f", 1);
  rv3_print(t, "t ", "%10.5f", 1);
  rv3_print(t, "t0", "%10.5f", 1);

  /* construct a configuration */
  for (i = 0; i < N; i++) {
    rv3_diff(tmp, y[i], t0);
    rm3_mulvec(x[i], r0, tmp);
    x[i][2] = -x[i][2];
  }
  
  printf("\n\nWith reflection\n");
  dev = rotfit3(x, NULL, y, NULL, N, r, t);  
  printf("dev = %g\n", dev);
  rm3_print(r, "R ", "%10.5f", 1);
  rm3_print(r0, "R0", "%10.5f", 1);
  rv3_print(t, "t ", "%10.5f", 1);
  rv3_print(t, "t0", "%10.5f", 1);

  return 0;
}
