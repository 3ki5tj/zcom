#include <string.h>
#include "lu.h"

#define N 3
int main(void)
{
  double A[N][N]={{1, 0.5, 1.0/3}, {0.5, 1.0/3, 0.25}, {1.0/3,.25,.2}}, a[N][N];
  double b[N]={1, 2, 3}, x[N];
  int i, j, index[N];

  memcpy(a, A, sizeof(double)*N*N);
  memcpy(x, b, sizeof(double)*N);
  lusolve((double *) a, x, N);
  for (i = 0; i < N; i++) 
    printf("%d: %g\n", i, x[i]);

  memcpy(a, A, sizeof(double)*N*N);
  memcpy(x, b, sizeof(double)*N);
  ludcmp((double *) a, index, N);
  for (i = 0; i < N; i++) {
    printf("%d: ", index[i]);
    for (j = 0; j < N; j++) {
      printf("%10.6f ", a[i][j]);
    }
    printf("\n");
  }
  lubksb((double *) a, x, index, N);
  for (i = 0; i < N; i++) 
    printf("%d: %g\n", i, x[i]);
  return 0;
}

