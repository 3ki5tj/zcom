#include "svd.c"

#define N 3

static void test_svdsolve(void)
{
  double a[N*N] = {1,2,3, 2,3,4, 2,4,5.9999}, b[N] = {2., 4, 4.0};
  double x[N], y;
  int i, j;

  svdsolve(a, x, b, N, 1e-6);

  printf("Equation:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%g ", a[i*N+j]);
    printf("\t%g\n", b[i]);
  }
  printf("x values:\n");
  for (i = 0; i < N; i++)
    printf("%g ", x[i]);
  printf("\n");
  for (i = 0; i < N; i++) {
    for (y = 0, j = 0; j < N; j++)
      y += a[i*N+j]*x[j];
    fprintf(stderr, "%g err = %g\n", b[i], y - b[i]);
  }
}

int main(void)
{
  test_svdsolve();
  return 0;
}
