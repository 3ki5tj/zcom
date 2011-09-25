#include <stdio.h>
#include <string.h>
#include "rv3.h"
#include "include/eig.c"

int main(void)
{
  real a[3][3] = {
    {2., 0, 0}, 
    {0, 1, 1},
    {0, 1, 1}}, mat[3][3];
  real v[3] = {0, 0, 0}, vecs[3][3];
  int j;
  
  /* cheap 3x3 eigensystem */ 
  memcpy(mat, a, sizeof(real)*9);
  rm3_eigval(v, mat);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  rm3_eigvecs(vecs, mat, v, 0);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
  
  /* general eigen system */
  printf("\n\nCHECKING...\n");
  memcpy(mat, a, sizeof(mat));
  eigsym((real *)mat, v, (real *)vecs, 3);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
  return 0;
}
