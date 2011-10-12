#include <stdio.h>
#include <string.h>
#define HAVEREAL 1
typedef float real;
#include "rv3.h"
#include "include/eig.c"

int main(void)
{
  real a[3][3] = {{2.f, 0.f, 0.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}}, mat[3][3];
  real v[3] = {0, 0, 0}, vecs[3][3];
  int i;

  for (i = 0; i < 50000000; i++) {
#ifdef CHEAP
    /* cheap 3x3 eigensystem */
    memcpy(mat, a, sizeof(real)*9);
    rm3_eigval(v, mat);
    rm3_eigvecs(vecs, mat, v, 0);
#else
    /* general eigen system */
    memcpy(mat, a, sizeof(mat));
    eigsym((real *)mat, v, (real *)vecs, 3);
#endif
  }
  return 0;
}
