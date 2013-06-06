#include "eig.h"

#define N  3

/* b = m * a */
void matmul(double *b, const double *m, const double *a, int n, int c)
{
  int i, j;
  for (i = 0; i < n; i++) {
    b[i] = 0.;
    for (j = 0; j < n; j++)
      b[i] += m[i*n + j]*a[j*n + c];
  }
}


int main(void)
{
  double mat[N*N], eig[N], vec[N*N], b[N];
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) mat[i*N+j] = 1.;
    mat[i*N+i] = i + 1.;
  }

  eigsym(mat, eig, vec, N);

  for (i = 0; i < N; i++) {
    printf("%d: %8.4f; ", i, eig[i]);
    /* print eigenvector i */
    for (j = 0; j < N; j++)
      printf("%8.4f ", vec[j*N+i]);
    printf("; ");
    /* verify eigenvectors */
    matmul(b, mat, vec, N, i);
    for (j = 0; j < N; j++)
      printf("%8.4f ", b[j]);
    printf("\n");
  }

  return 0;
}
