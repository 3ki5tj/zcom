
#define HAVEREAL 1
typedef float real;

#include "rv3.c"
#include "eig.h"
#include <time.h>

int simple = 1;   /* simple test of a single matrix a */
int sptest = 1;   /* random matrix speed test */
int rndtest = 1;  /* random matrix test */


real a[3][3] = {{3.f, 0.f, 0.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}};


static void rm3_rndsym0(real m[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) {
      m[i][j] = (real) rnd0();
      if (j != i) m[j][i] = m[i][j];
    }
}


static void testsimp(void)
{
  real v[3] = {0, 0, 0}, vecs[3][3], mat[3][3];
  int j;

  /* cheap 3x3 eigensystem */
  memcpy(mat, a, sizeof(real)*9);
  rm3_eigsys(v, vecs, mat, 0);
  printf("eigenvalues:   %20.12f, %20.12f, %20.12f\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %20.12f, %20.12f, %20.12f\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);

  /* general eigen system */
  printf("\n\nCHECKING using eigsys()...\n");
  memcpy(mat, a, sizeof(mat));
  eigsym((real *) mat, v, (real *) vecs, 3);
  printf("eigenvalues:   %20.12f, %20.12f, %20.12f\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %20.12f, %20.12f, %20.12f\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
}


static void testspeed(int ntests)
{
  real mat[3][3], v[3] = {0, 0, 0}, vecs[3][3];
  int i;
  clock_t t0;

  t0 = clock();
  for (i = 0; i < ntests; i++) {
    /* cheap 3x3 eigensystem */
    rm3_rndsym0(mat);
    rm3_eigsys(v, vecs, mat, 0);
  }
  printf("rm3_eigsys() used %fs\n", (double) (clock() - t0) / CLOCKS_PER_SEC);

  t0 = clock();
  for (i = 0; i < ntests; i++) {
    /* general eigen system */
    rm3_rndsym0(mat);
    eigsym((real *) mat, v, (real *) vecs, 3);
  }
  printf("eigsys() used %fs\n", (double) (clock() - t0) / CLOCKS_PER_SEC);
}


/* random symmetrical matrix */
static void testrnd(int nrands)
{
  int t;
  real maxdet = 0, maxnorm = 0, maxeig = 0;

  const double tol = 1e-3;

  real m[3][3], vals[3], vecs[3][3], vs[3][3], det, del, vn;
  int i;

  for (t = 0; t < nrands; t++) {
    /* construct a random symmetric matrix */
    rm3_rndsym0(m);

    /* compute eigenvalues and eigenvectors */
    rm3_eigsys(vals, vecs, m, 1);

    /* compute the determinant of the matrix of eigenvectors */
    det = rm3_det(vecs);
    del = (real) fabs(fabs(det) - 1);
    if (del > maxdet) maxdet = del;
    if (del > tol) {
      fprintf(stderr, "fatal: det %g\n", det);
      goto ERR;
    }

    /* accuracy of eigenvalues */
    for (i = 0; i < 3; i++) {
      rm3_mulvec(vs[i], m, vecs[i]);

      /* check the magnitude of vs */
      vn = rv3_norm(vs[i]);
      del = (real) fabs(vn - fabs(vals[i]));
      if (del > maxnorm) maxnorm = del;
      if (del > tol) {
        fprintf(stderr, "fatal: eigenvalue %d: norm(vs) = %g vs %g\n", i, vn, vals[i]);
        goto ERR;
      }

      /* check if vs[i] is parallel to vecs[i] */
      vn = rv3_dot(vs[i], vecs[i]);
      del = (real) fabs(vn - vals[i]);
      if (del > maxeig) maxeig = del;
      if (del > tol) {
        fprintf(stderr, "fatal: eigenvalue %d: vs.vec %g vs. %g\n", i, vn, vals[i]);
        goto ERR;
      }
    }
  }
  printf("pass %d tests: maxdet %g, maxnorm %g, maxeig %g\n", nrands, maxdet, maxnorm, maxeig);
  return;
ERR:
  rm3_print(m,    "m",    "%20.12f", 1);
  rm3_print(vecs, "vecs", "%20.12f", 1);
  rv3_print(vals, "vals", "%20.12f", 1);
  return;
}

int main(void)
{
  if (simple) testsimp();
  if (sptest) testspeed(2000000);
  if (rndtest) testrnd(1000000);
  return 0;
}

