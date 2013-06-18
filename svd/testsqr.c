#include "svd.h"

#define N 3

static void test_svd(void)
{
  double a[N*N] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double u[N*N], v[N*N], w[N];
  int i, j;

  for (i = 0; i < N*N; i++) u[i] = a[i];
  if (svd(u, w, v, N, N) != 0) {
    fprintf(stderr, "svd failed\n");
    return;
  }
  printf("U matrix:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%g ", u[i*N+j]);
    printf("\n");
  }
  printf("S values:\n");
  for (i = 0; i < N; i++)
    printf("%g ", w[i]);
  printf("\n");
  printf("V matrix:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%g ", v[i*N+j]);
    printf("\n");
  }
}

int main(void)
{
  test_svd();
  return 0;
}
