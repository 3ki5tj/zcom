#include <stdio.h>
#include <string.h>
/*
#define HAVEREAL 1
typedef float real;
*/

#define RV3_DEBUG

#include "rv3.h"
#include "include/eig.c"

/*
real a[3][3] = {{2.f, 0.f, 0.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}};
*/
/*
real a[3][3] = {{17.f, 22.f, 0.f}, {22.f, 29.f, 0.f}, {0.f, 0.f, 81.f}};
*/

real a[3][3] = {
{  1.8771876912942007e-07, -6.5304014581726230e-06, -9.6244241206405254e-05},
{ -6.5304014581726230e-06,  2.2718102938071863e-04,  3.3481656417729158e-03},
{ -9.6244241206405254e-05,  3.3481656417729158e-03,  4.9344847126131020e-02}};

int main(void)
{
  real v[3] = {0, 0, 0}, vecs[3][3], mat[3][3];
  int j;
  
  /* cheap 3x3 eigensystem */ 
  memcpy(mat, a, sizeof(real)*9);
  rm3_eigsys(v, vecs, mat, 0);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
  
  /* general eigen system */
  printf("\n\nCHECKING using eigsys()...\n");
  memcpy(mat, a, sizeof(mat));
  eigsym((real *)mat, v, (real *)vecs, 3);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
  return 0;
}
