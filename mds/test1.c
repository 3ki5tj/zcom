#include "mds.c"

double *mkdismat0(void)
{
  double *dm;
  int i, n = 3;

  xnew(dm, n*n);
  for (i = 0; i < n; i++) dm[i*n+i] = 0.;
  dm[0*n+1] = dm[1*n+0] = 2.;
  dm[2*n+1] = dm[1*n+2] = 2.;
  dm[0*n+2] = dm[2*n+0] = 2.;
  return dm;
}

int main(void)
{
  int dim = 2;
  int n = 3;
  int i, d;
  double *dm, *x;

  xnew(x, n*dim); /* coordinates */
  dm = mkdismat0();
  mds_min0(x, dm, n, dim, 1e-14);
  for (i = 0; i < n; i++) {
    printf("%d: ", i);
    for (d = 0; d < dim; d++)
      printf("%g ", x[i*dim+d]);
    printf("\n");
  }
  free(x);
  return 0;
}
