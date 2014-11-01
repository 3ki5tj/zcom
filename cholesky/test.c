#include <string.h>
#include "cholesky.h"

#define N 3
int main(void)
{
  real a[N][N] = {{1, 2, 1}, {2, 5, 4}, {1, 4, 6}}, b[N] = {1, 2, 3};
  real l[N][N], x[N], y[N], z[N][N], az[N][N];
  int i, j, k;

  memcpy(l, a, N*N*sizeof(a[0][0]));
  memcpy(x, b, N*sizeof(b[0]));
  cholsolve((real *) l, x, N);

  printf("                  A                         b\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%10.6f  ", a[i][j]);
    printf(" |  %10.6f \n", b[i]);
  }

  printf("                  L\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%10.6f  ", l[i][j]);
    printf("\n");
  }

  printf("Solution x:\n");
  for ( i = 0; i < N; i++ )
    printf("%10.6f  ", x[i]);
  printf("\n");

  for ( i = 0; i < N; i++ )
    for ( y[i] = 0, j = 0; j < N; j++ )
      y[i] += a[i][j] * x[j];

  printf("      A x        b\n");
  for ( i = 0; i < N; i++ )
    printf("%10.6f  %10.6f\n", y[i], b[i]);
  printf("\n");

  memcpy(l, a, N*N*sizeof(a[0][0]));
  cholinv((real *) l, (real *) z, N);
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
