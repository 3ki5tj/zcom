#include <string.h>
#include "lu.h"

#define N 3
int main(void)
{
  double a[N][N]={{1, 0.5, 1.0/3}, {0.5, 1.0/3, 0.25}, {1.0/3,.25,.2}}, lu[N][N];
  double b[N]={1, 2, 3}, x[N], y[N], z[N][N], az[N][N];
  int i, j, k;

  memcpy(lu, a, sizeof(double)*N*N);
  memcpy(x, b, sizeof(double)*N);
  lusolve((double *) lu, x, N);
  printf("equation:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%10.6f ", a[i][j]);
    printf("\t%10.6f\n", b[i]);
  }
  printf("matrix LU:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%10.6f ", lu[i][j]);
    printf("\n");
  }
  printf("solution x:\n");
  for (i = 0; i < N; i++)
    printf("%d: %g\n", i, x[i]);

  for (i = 0; i < N; i++)
    for (y[i] = 0, j = 0; j < N; j++)
      y[i] += a[i][j]*x[j];
  printf("        A x        b\n");
  for (i = 0; i < N; i++)
    printf("%d: %10.6f %10.6f\n", i, y[i], b[i]);

  memcpy(lu, a, N*N*sizeof(a[0][0]));
  luinv((real *) lu, (real *) z, N);
  printf("                A^(-1)\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%10.6f  ", z[i][j]);
    printf("\n");
  }

  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ ) {
      az[i][j] = 0;
      for ( k = 0; k < N; k++ )
        az[i][j] += a[i][k] * z[k][j];
    }
  }
  printf("                A.A^(-1)\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%10.6f  ", az[i][j]);
    printf("\n");
  }

  return 0;
}

