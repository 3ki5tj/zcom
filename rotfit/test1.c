#include "rotfit.c"

#define N 5
int main(void)
{
  int i, j;
  real x[N][3], y[N][3], tmp[3], dev;
  real r0[3][3] = {{1, 0, 0}, {0, 0.866, 0.5}, {0, -.5, .866}}, 
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
    mat3_mulvec(x[i], r0, tmp);
  }
  
  printf("\n\nPure rotation\n");
  rotfit(x, y, NULL, N, r, t);
  dev = rmsdev(x, NULL, y, NULL, N, r, t);  
  printf("dev = %g\n", dev);
  printf("r  =\n %g %g %g\n %g %g %g\n %g %g %g\n", r[0][0],r[0][1],r[0][2],
      r[1][0],r[1][1],r[1][2],r[2][0],r[2][1],r[2][2]);
  printf("r0 =\n %g %g %g\n %g %g %g\n %g %g %g\n", r0[0][0],r0[0][1],r0[0][2],
      r0[1][0],r0[1][1],r0[1][2],r0[2][0],r0[2][1],r0[2][2]);
  printf("t = %g %g %g, t0 = %g %g %g\n", t[0],t[1],t[2],t0[0],t0[1],t0[2]);

  /* construct a configuration */
  for (i = 0; i < N; i++) {
    rv3_diff(tmp, y[i], t0);
    mat3_mulvec(x[i], r0, tmp);
    x[i][2] = -x[i][2];
  }
  
  printf("\n\nWith reflection\n");
  rotfit(x, y, NULL, N, r, t);
  dev = rmsdev(x, NULL, y, NULL, N, r, t);  
  printf("dev = %g\n", dev);
  printf("r  =\n %g %g %g\n %g %g %g\n %g %g %g\n", r[0][0],r[0][1],r[0][2],
      r[1][0],r[1][1],r[1][2],r[2][0],r[2][1],r[2][2]);
  printf("r0 =\n %g %g %g\n %g %g %g\n %g %g %g\n", r0[0][0],r0[0][1],r0[0][2],
      r0[1][0],r0[1][1],r0[1][2],r0[2][0],r0[2][1],r0[2][2]);
  printf("t = %g %g %g, t0 = %g %g %g\n", t[0],t[1],t[2],t0[0],t0[1],t0[2]);

  return 0;
}
