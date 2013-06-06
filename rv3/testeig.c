#include <stdio.h>
#include <string.h>

#define HAVEREAL 1
typedef float real;

#include "rv3.h"
#include "include/eig.h"

int simple = 1;   /* simple test of a single matrix a */
int sptest = 0;   /* do speed test */
int rndtest = 0;  /* random matrix test */

/*
real a[3][3] = {{2.f, 0.f, 0.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}};
*/
/*
real a[3][3] = {{17.f, 22.f, 0.f}, {22.f, 29.f, 0.f}, {0.f, 0.f, 81.f}};
*/

real a[3][3] = {
{  1.8771876912942007e-07, -6.5304014581726230e-06, -9.6244241206405254e-05},
{ -6.5304014581726230e-06,  2.2718102938071863e-04,  3.3481656417729158e-03},
{ -9.6244241206405254e-05,  3.3481656417729158e-03,  4.9344847126131020e-02}};

static void testsimp(void)
{
  real v[3] = {0, 0, 0}, vecs[3][3], mat[3][3];
  int j;

  /* cheap 3x3 eigensystem */
  memcpy(mat, a, sizeof(real)*9);
  rm3_eigsys(v, vecs, mat, 0);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);

  /* general eigen system */
  printf("\n\nCHECKING using eigsys()...\n");
  memcpy(mat, a, sizeof(mat));
  eigsym((real *) mat, v, (real *) vecs, 3);
  printf("eigenvalues are %g, %g, %g\n", v[0], v[1], v[2]);
  for (j = 0; j < 3; j++)
    printf("eigenvector %d: %g, %g, %g\n", j, vecs[0][j], vecs[1][j], vecs[2][j]);
}

static void testspeed(int ntests)
{
  real a[3][3] = {{2.f, 0.f, 0.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}}, mat[3][3];
  real v[3] = {0, 0, 0}, vecs[3][3];
  int i;

  for (i = 0; i < ntests; i++) {
#ifdef CHEAP
    /* cheap 3x3 eigensystem */
    memcpy(mat, a, sizeof(real)*9);
    rm3_eigsys(v, vecs, mat, 0);
#else
    /* general eigen system */
    memcpy(mat, a, sizeof(mat));
    eigsym((real *) mat, v, (real *) vecs, 3);
#endif
  }
}

/* random symmetrical matrix */
static void testrnd(int nrands)
{
  int t;
  real maxdet = 0, maxnorm = 0, maxeig = 0;

  const double tol = 1e-4;

  real m[3][3], vals[3], vecs[3][3], vs[3][3], det;
  real del;
  int i, j;

  for (t = 0; t < nrands; t++) {
    /* construct a random symmetric matrix */
    for (i = 0; i < 3; i++)
      for (j = i; j < 3; j++) {
        m[i][j] = 1.0 * rand()/RAND_MAX;
        if (j != i) m[j][i] = m[i][j];
      }

    /* compute eigenvalues and eigenvectors */
    rm3_eigsys(vals, vecs, m, 1);

    /* compute det */
    det = rm3_det(vecs);
    del = fabs(fabs(det) - 1);
    if (del > maxdet) maxdet = del;
    if (del > tol) {
      fprintf(stderr, "fatal: det %g\n", det);
      goto ERR;
    }

    /* accuracy of eigenvalues */
    for (i = 0; i < 3; i++) {
      real vn;

      rm3_mulvec(vs[i], m, vecs[i]);

      /* check the magnitude of vs */
      vn = rv3_norm(vs[i]);
      del = fabs(vn - fabs(vals[i]));
      if (del > maxnorm) maxnorm = del;
      if (del > tol) {
        fprintf(stderr, "fatal; eigenvalue %d: norm(vs) = %g vs %g\n", i, vn, vals[i]);
        rv3_print(vecs[i],  "vecs[i]",  "%24.14e", 1);
        rv3_print(vs[i],    "vs[i]",    "%24.14e", 1);
        goto ERR;
      }

      /* check if vs[i] is parallel to vecs[i] */
      vn = rv3_dot(vs[i], vecs[i]);
      del = fabs(vn - vals[i]);
      if (del > maxeig) maxeig = del;
      if (del > tol) {
        fprintf(stderr, "fatal: eigenvalue %d: vs.vec %g vs. %g\n", i, vn, vals[i]);
        rv3_print(vecs[i],  "vecs[i]",  "%24.14e", 1);
        rv3_print(vs[i],    "vs[i]",    "%24.14e", 1);
        goto ERR;
      }
    }
  }
  printf("pass %d tests: maxdet %g, maxnorm %g, maxeig %g\n", nrands, maxdet, maxnorm, maxeig);
  return;
ERR:
  rm3_print(m,    "m",    "%24.14e", 1);
  rm3_print(vecs, "vecs", "%24.14e", 1);
  rv3_print(vals, "vals", "%24.14e", 1);
  exit(1);
}

int main(void)
{
  if (simple) testsimp();
  if (sptest) testspeed(1000000);
  if (rndtest) testrnd(1000000);
  return 0;
}

