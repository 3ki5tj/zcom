#include "rv3.h"
#ifndef ROTFIT_C__
#define ROTFIT_C__
#include "rotfit.h"

/* least square fit x to y after rotation/translation */
real rotfit3(rv3_t *x, rv3_t *xf, rv3_t *y, const real *w, int n,
    real (*r)[3], real *t)
{
  int i;
  real wtot = 0, sq, dev = 0;
  rv3_t xc, yc, xs, ys, sig, t_;
  real u[3][3], v[3][3], s[3][3] = {{0,0,0},{0,0,0},{0,0,0}}, xy[3][3], r_[3][3];

  if (r == NULL) r = r_;
  if (t == NULL) t = t_;
    
  /* 1. compute the centers */
  rv3_zero(xc);
  rv3_zero(yc);
  if (w == NULL) {
    for (i = 0; i < n; i++) {
      rv3_inc(xc, x[i]);
      rv3_inc(yc, y[i]);
    }
    wtot = n;
  } else {
    for (wtot = 0., i = 0; i < n; i++) {
      rv3_sinc(xc, x[i], w[i]);
      rv3_sinc(yc, y[i], w[i]);
      wtot += w[i];
    }
  }
  rv3_smul(xc, 1.f/wtot);
  rv3_smul(yc, 1.f/wtot);

  /* 2. compute 3x3 asymmetric covarience matrix S = (x-xc) (y-yc)^T */
  for (i = 0; i < n; i++) {
    rv3_diff(xs, x[i], xc); /* shift to the center avoid the translation */
    rv3_diff(ys, y[i], yc);
    mat3_vtv(xy, xs, ys);
    sq  = rv3_sqr(xs);
    sq += rv3_sqr(ys);
    if (w) {
      mat3_sinc(s, xy, w[i]);
      dev += w[i]*sq;
    } else {      
      mat3_inc(s, xy);
      dev += sq; /* Tr(x^T x + y^T y) */
    }
  }

  /* 3. SVD decompose S = u sig v^T */
  mat3_svd(s, u, sig, v);

  /* 4. compute R = v u^T */
  mat3_mult(r, v, u);
  if (mat3_det(r) < 0) { /* to avoid a reflection */
    mat3_trans(u);
    rv3_neg(u[2]); /* flip the last eigenvalue */
    mat3_mul(r, v, u);
    dev -= 2*(sig[0]+sig[1]-sig[2]);
  } else {
    dev -= 2*(sig[0]+sig[1]+sig[2]); /* -2 Tr(R x y^T) */
  }
  rv3_diff(t, yc, mat3_mulvec(xs, r, xc)); /* t = yc - R xc */

  /* 5. compute the rotated structure */
  if (xf) {
    for (dev = 0, i = 0; i < n; i++) {
      rv3_add(xf[i], mat3_mulvec(xs, r, x[i]), t); /* xf = R x + t */
      sq = rv3_dist2(y[i], xf[i]);
      dev +=  (w ? w[i]*sq : sq); /* recompute the deviation */
    }
  }
  return (real) sqrt(dev/wtot);
}

#endif

