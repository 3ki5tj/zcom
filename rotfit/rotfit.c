#include "util.h"
#include "rv3.h"
#ifndef ROTFIT_C__
#define ROTFIT_C__
#include "rotfit.h"

/* rotate x to fit into y, as r x + t */
static int rotfit3_rt(rv3_t *x0, rv3_t *y0, const real *w, int n,
    real r[3][3], real t[3])
{
  int i, id, jd;
  real wtot = 0.f, sm, tmp;
  rv3_t *x, *y, xc, yc, sig;
  real u[3][3], v[3][3], s[3][3];

  /* 1. compute the centers */
  rv3_zero(xc);
  rv3_zero(yc);
  if (w == NULL) {
    for (i = 0; i < n; i++) {
      rv3_inc(xc, x0[i]);
      rv3_inc(yc, y0[i]);
    }
    wtot = n;
  } else {
    for (wtot = 0., i = 0; i < n; i++) {
      rv3_sinc(xc, x0[i], w[i]);
      rv3_sinc(yc, y0[i], w[i]);
      wtot += w[i];
    }
  }
  wtot = 1./wtot;
  rv3_smul(xc, wtot);
  rv3_smul(yc, wtot);
  
  /* 2. compute centered vectors */
  xnew(x, n);
  xnew(y, n);
  for (i = 0; i < n; i++) {
    rv3_diff(x[i], x0[i], xc);
    rv3_diff(y[i], y0[i], yc);
  }

  /* 3. compute 3x3 covarience matrix */
  for (id = 0; id < 3; id++)
    for (jd = 0; jd < 3; jd++) {
      for(sm = 0, i = 0; i < n; i++) {
        tmp = x[i][id]*y[i][jd];
        if (w) sm += w[i]*tmp;
        else sm += tmp;
      }
      s[id][jd] = sm;
    }

  /* 4. SVD decompose S, compute v u^T 
   * note s is not symmetric, SVD is necessary */
  mat3_svd(s, u, sig, v);
  /* R = v u^T */ 
  mat3_mult(r, v, u);
  if (mat3_det(r) < 0) { /* reflection */
    mat3_trans(u);
    rv3_neg(u[2]);
    mat3_mul(r, v, u);
  }
  rv3_copy(t, yc);
  mat3_mulvec(yc, r, xc);
  rv3_dec(t, yc);
  free(x);
  free(y);
  return 0;
}

/* calculate RMS difference */
static real rotfit3_rmsd(rv3_t *x0, rv3_t *x, rv3_t *y0, const real *w, int n,
    real r[3][3], const real t[3])
{
  int i, hasx = 1;
  real dev = 0, sq, wtot = 0, rx[3];

  if (x == NULL) {
    hasx = 0;
    xnew(x, n);
  }
  for (i = 0; i < n; i++) {
    mat3_mulvec(rx, r, x0[i]);
    rv3_add(x[i], rx, t);
    sq = rv3_dist2(y0[i], x[i]);
    if (w) {
      dev += w[i]*sq;
      wtot += w[i];
    } else {
      dev += sq;
      wtot += 1.;
    }
  }
  if (!hasx) free(x);
  return sqrt(dev/n);
}

/* fit x0 to y after rotation/translation */
real rotfit3(rv3_t *x0, rv3_t *x, rv3_t *y0, const real *w, int n, 
    real (*r)[3], real *t)
{
  unsigned flags = 0;
  real rmsd;
  if (r == NULL) { xnew(r, 3); flags |= 1; }
  if (t == NULL) { xnew(t, 3); flags |= 2; }
  rotfit3_rt(x0, y0, w, n, r, t);
  rmsd = rotfit3_rmsd(x0, x, y0, w, n, r, t);
  if (flags & 1) free(r);
  if (flags & 2) free(t);
  return rmsd;
}

#endif

