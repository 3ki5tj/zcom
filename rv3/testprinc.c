#include "rv3.c"

#define N 20



/* make a helcial polypeptide */
static void mkhelix(real x[][3], int n)
{
  double ang = 0., z = 0.;
  int i;

  for (i = 0; i < n; i++) {
    x[i][0] = cos(ang);
    x[i][1] = sin(ang);
    x[i][2] = z;
    ang += 100./180.*M_PI;
    z += 1.25;
    printf("%2d: %8.3f %8.3f %8.3f\n", i, x[i][0], x[i][1], x[i][2]);
  }
}


static void calcmmt(real mat[3][3], real xc[3], real x[][3], int n)
{
  real dx[3];
  int i, j, k;

  rv3_zero(xc);
  for (i = 0; i < n; i++)
    rv3_inc(xc, x[i]);
  rv3_smul(xc, 1./n);
  for (j = 0; j < 3; j++)
    for (k = 0; k < 3; k++)
      mat[j][k] = 0.;
  for (i = 0; i < n; i++) {
    rv3_diff(dx, x[i], xc);
    for (j = 0; j < 3; j++)
      for (k = j; k < 3; k++)
        mat[j][k] += dx[j]*dx[k];
  }
  for (j = 0; j < 3; j++) {
    for (k = j+1; k < 3; k++)
      mat[k][j] = mat[j][k] = mat[j][k]/n;
    mat[j][j] *= 1./n;
  }
}


int main(void)
{
  real x[N][3] = {{0., 0., 0.}}, mat[3][3], xc[3];
  real v[3] = {0, 0, 0}, vecs[3][3] = {{0}};
  real x0[3], x1[3];

  mkhelix(x, N);
  calcmmt(mat, xc, x, N);
  rm3_eigval(v, mat);
  rm3_eigvecs(vecs, mat, v[0]);
  printf("val: %g, %g, %g; vec = {%g, %g, %g}\n",
      v[0], v[1], v[2], vecs[0][0], vecs[0][1], vecs[0][2]);
  printf("xc = %g, %g, %g\n", xc[0], xc[1], xc[2]);
  rv3_sadd(x0, xc, vecs[0], -sqrt(v[0]) * sqrt(3));
  rv3_sadd(x1, xc, vecs[0], +sqrt(v[0]) * sqrt(3));
  printf("%g, %g, %g;   %g, %g, %g\n", x0[0], x0[1], x0[2], x1[0], x1[1], x1[2]);
  return 0;
}

