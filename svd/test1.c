#include "svd.c"

#define M 3
#define N 2

static void test_svd(void)
{
  double a[M*N] = {1, 2,  2, 3,  3, 4};
  double u[M*N], v[N*N], w[N];
  int i, j;

  for (i = 0; i < M*N; i++) u[i] = a[i];
  if (svd(u, w, v, M, N) != 0) {
    fprintf(stderr, "svd failed\n");
    return;
  }
  printf("U matrix:\n");
  for (i = 0; i < M; i++) {
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
