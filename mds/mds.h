#include "util.h"
#include "eig.h"
#ifndef MDS_H__
#define MDS_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* the main function for multidimensional scaling */
real mds_min0(real *x, real *dm, int n, int dim, double tol);

/* compute force and energy */
static real mds_force(real *x, real *f, real *dm, int n, int dim)
{
  const real dmin = 1e-6f;
  int i, j, k;
  real ene = 0., dij, dref, dsc;
  real *dx, *xi, *xj, *fi, *fj;

  xnew(dx, dim);
  for (i = 0; i < n*dim; i++) f[i] = 0.;
  for (i = 0; i < n; i++) {
    xi = x + i*dim;
    fi = f + i*dim;
    for (j = i+1; j < n; j++) {
      xj = x + j*dim;
      fj = f + j*dim;
      for (dij = 0, k = 0; k < dim; k++) {
        dx[k] = xi[k] - xj[k];
        dij += dx[k]*dx[k];
      }
      dij = (real) sqrt(dij);
      dref = dm[i*n+j];
      if (dref < dmin) dref = dmin;
      dsc = dij/dref - 1;
      ene += dsc*dsc;
      /* dij is to be used in denominator in the loop, ensure its > 0 */
      if (dij < dmin) dij = dmin;
      for (k = 0; k < dim; k++) {
        dx[k] *= -(dij - dref)/(dref*dref*dij);
        fi[k] += dx[k];
        fj[k] -= dx[k];
      }
    }
  }
  free(dx);
  return ene;
}

/* make coordinates neat
 * center coordinates
 * rotate to principle coordinates */
static void mds_trim(real *x, int n, int dim)
{
  real *av, *mmt, *eig, *vec, *xi, *b;
  int i, d, d2;

  /* center the graph */
  xnew(av, dim);
  for (i = 0; i < n; i++)
    for (d = 0; d < dim; d++)
      av[d] += x[i*dim + d];
  for (d = 0; d < dim; d++) av[d] /= n;
  for (i = 0; i < n; i++)
    for (d = 0; d < dim; d++)
      x[i*dim + d] -= av[d];
  free(av);

  /* compute moments */
  xnew(mmt, dim*dim);
  xnew(eig, dim);
  xnew(vec, dim*dim);
  xnew(b, dim);
  for (i = 0; i < n; i++) {
    for (d = 0; d < dim; d++)
      for (d2 = 0; d2 < dim; d2++)
        mmt[d*dim+d2] += x[i*dim+d]*x[i*dim+d2];
  }
  eigsym(mmt, eig, vec, dim);
  for (i = 0; i < n; i++) {
    /* rotate x[i*dim .. i*dim+d-1] --> b[0 .. d-1] */
    xi = x + i*dim;
    for (d = 0; d < dim; d++) /* b = xi.vec */
      for (b[d] = 0., d2 = 0; d2 < dim; d2++)
        b[d] += xi[d2]*vec[d2*dim+d];
    for (d = 0; d < dim; d++) xi[d] = b[d];
  }
  free(b);
  free(eig);
  free(vec);
  free(mmt);
}

/* multidimensional scaling - steepest descend
 * given a distance matrix dm[n x n],
 * return best mds position x[n x dim];
 * dim is the target dimensional, e.g. 2
 * return the total discrepancy */
real mds_min0(real *x, real *dm, int n, int dim, double tol)
{
  int i, j, it, itermax = 100000, npr;
  real *f, *fp, *xp, ene, enep;
  real dt = 1e-1f;

  if (n == 1) {
    for (j = 0; j < dim; j++) x[j] = 0.;
    return 0.0;
  }
  npr = n*(n-1)/2;
  xnew(f, n*dim);
  xnew(xp, n*dim);
  xnew(fp, n*dim);
  for (i = 0; i < n*dim; i++)
    x[i] = 1.f*rand()/RAND_MAX;
  ene = mds_force(x, f, dm, n, dim);
  for (it = 0; it < itermax; it++) {
    enep = ene;
    for (i = 0; i < n*dim; i++) { /* backup */
      xp[i] = x[i];
      fp[i] = f[i];
    }
    for (i = 0; i < n*dim; i++) x[i] += f[i]*dt;
    ene = mds_force(x, f, dm, n, dim);
    if (ene > enep) {
      dt *= 0.7f;
      for (i = 0; i < n*dim; i++) { /* recover */
        x[i] = xp[i];
        f[i] = fp[i];
      }
    } else {
      if (fabs(ene-enep) < tol*npr*dt) break;
      dt *= 1.03f; /* attempt to increase time step */
    }
  }
  if (it >= itermax) {
    fprintf(stderr, "mds: failed to converge after %d iterations, %g\n",
        it, fabs(ene-enep));
  }
  mds_trim(x, n, dim);
  free(xp);
  free(f);
  free(fp);
  return ene;
}

#endif /* MDS_H__ */

