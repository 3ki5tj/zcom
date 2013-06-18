#include "rv3.h"
#ifndef RM3_H__
#define RM3_H__
/* routines for 3x3 matrices */


INLINE dv3_t *dm3_fromrm3(double a[3][3], real ra[3][3])
{
  a[0][0] = ra[0][0];
  a[0][1] = ra[0][1];
  a[0][2] = ra[0][2];
  a[1][0] = ra[1][0];
  a[1][1] = ra[1][1];
  a[1][2] = ra[1][2];
  a[2][0] = ra[2][0];
  a[2][1] = ra[2][1];
  a[2][2] = ra[2][2];
  return a;
}


INLINE rv3_t *rm3_fromdm3(real a[3][3], double da[3][3])
{
  a[0][0] = (real) da[0][0];
  a[0][1] = (real) da[0][1];
  a[0][2] = (real) da[0][2];
  a[1][0] = (real) da[1][0];
  a[1][1] = (real) da[1][1];
  a[1][2] = (real) da[1][2];
  a[2][0] = (real) da[2][0];
  a[2][1] = (real) da[2][1];
  a[2][2] = (real) da[2][2];
  return a;
}



/* a = b */
INLINE rv3_t *rm3_copy(real a[3][3], real b[3][3])
{
  rv3_copy(a[0], b[0]);
  rv3_copy(a[1], b[1]);
  rv3_copy(a[2], b[2]);
  return a;
}



/* transpose */
INLINE rv3_t *rm3_trans(real a[3][3])
{
  real x;

  x = a[0][1], a[0][1] = a[1][0], a[1][0] = x;
  x = a[0][2], a[0][2] = a[2][0], a[2][0] = x;
  x = a[2][1], a[2][1] = a[1][2], a[1][2] = x;
  return a;
}



/* a = u^T v */
INLINE rv3_t *rm3_vtv(real a[3][3], const real *u, const real *v)
{
  a[0][0] = u[0] * v[0];
  a[0][1] = u[0] * v[1];
  a[0][2] = u[0] * v[2];
  a[1][0] = u[1] * v[0];
  a[1][1] = u[1] * v[1];
  a[1][2] = u[1] * v[2];
  a[2][0] = u[2] * v[0];
  a[2][1] = u[2] * v[1];
  a[2][2] = u[2] * v[2];
  return a;
}



/* a += b */
INLINE rv3_t *rm3_inc(real a[3][3], real b[3][3])
{
  a[0][0] += b[0][0];
  a[0][1] += b[0][1];
  a[0][2] += b[0][2];
  a[1][0] += b[1][0];
  a[1][1] += b[1][1];
  a[1][2] += b[1][2];
  a[2][0] += b[2][0];
  a[2][1] += b[2][1];
  a[2][2] += b[2][2];
  return a;
}



/* a += b*s */
INLINE rv3_t *rm3_sinc(real a[3][3], real b[3][3], real s)
{
  a[0][0] += b[0][0] * s;
  a[0][1] += b[0][1] * s;
  a[0][2] += b[0][2] * s;
  a[1][0] += b[1][0] * s;
  a[1][1] += b[1][1] * s;
  a[1][2] += b[1][2] * s;
  a[2][0] += b[2][0] * s;
  a[2][1] += b[2][1] * s;
  a[2][2] += b[2][2] * s;
  return a;
}



/* c = a b */
INLINE rv3_t *rm3_mul(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
  return c;
}



/* c = a b^T */
INLINE rv3_t *rm3_mult(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = rv3_dot(a[i], b[j]);
  return c;
}



/* c = a^T b */
INLINE rv3_t *rm3_tmul(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[0][i]*b[0][j] + a[1][i]*b[1][j] + a[2][i]*b[2][j];
  return c;
}



/* c = a v */
INLINE real *rm3_mulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0] + a[0][1]*v[1] + a[0][2]*v[2];
  c[1] = a[1][0]*v[0] + a[1][1]*v[1] + a[1][2]*v[2];
  c[2] = a[2][0]*v[0] + a[2][1]*v[1] + a[2][2]*v[2];
  return c;
}



/* c = a^T v */
INLINE real *rm3_tmulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0] + a[1][0]*v[1] + a[2][0]*v[2];
  c[1] = a[0][1]*v[0] + a[1][1]*v[1] + a[2][1]*v[2];
  c[2] = a[0][2]*v[0] + a[1][2]*v[1] + a[2][2]*v[2];
  return c;
}



/* return 0 rotation matrix around v for ang */
INLINE rv3_t *rm3_mkrot(real m[3][3], const real *v, real ang)
{
  real c = (real) cos(ang), s = (real) sin(ang), nc, n[3];

  rv3_copy(n, v);
  rv3_normalize(n);
  nc = 1 - c;
  m[0][0] = n[0] * n[0] * nc + c;
  m[0][1] = n[0] * n[1] * nc - n[2] * s;
  m[0][2] = n[0] * n[2] * nc + n[1] * s;
  m[1][0] = n[1] * n[0] * nc + n[2] * s;
  m[1][1] = n[1] * n[1] * nc + c;
  m[1][2] = n[1] * n[2] * nc - n[0] * s;
  m[2][0] = n[2] * n[0] * nc - n[1] * s;
  m[2][1] = n[2] * n[1] * nc + n[0] * s;
  m[2][2] = n[2] * n[2] * nc + c;
  return m;
}



/* rotate v0 around u by ang, save result to v1 */
INLINE real *rv3_rot(real *v1, const real *v0, const real *u, real ang)
{
  real m[3][3];

  rm3_mkrot(m, u, ang);
  rm3_mulvec(v1, m, v0);
  return v1;
}



/* determinant of a 3x3 matrix */
INLINE real rm3_det(real a[3][3])
{
  return a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
      +  a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a[2][2])
      +  a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]);
}



/* inverse matrix b = a^(-1) */
INLINE rv3_t *rm3_inv(real b[3][3], real a[3][3])
{
  real d00, d01, d02, detm;
  d00 = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  d01 = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  d02 = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  detm = a[0][0]*d00 + a[0][1]*d01 + a[0][2]*d02;
  if (fabs(detm) < 1e-30) detm = (detm < 0) ? -1e-30f: 1e-30f;
  b[0][0] = d00/detm;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detm;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/detm;
  b[1][0] = d01/detm;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/detm;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/detm;
  b[2][0] = d02/detm;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/detm;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detm;
  return b;
}



/* compute eigenvalues of a 3x3 matrix
 * solving a cubic equation
 * use double for internal calculation */
INLINE real *rm3_eigval(real v[3], real a[3][3])
{
  double m, p, q, pr, pr3, a00, a11, a22;

  m = (a[0][0] + a[1][1] + a[2][2])/3;
  a00 = a[0][0] - m;
  a11 = a[1][1] - m;
  a22 = a[2][2] - m;
  q = ( a00 * (a11*a22 - a[1][2]*a[2][1])
      + a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a22)
      + a[0][2] * (a[1][0]*a[2][1] - a11*a[2][0]) ) / 2.0;
  p = (a00*a00 + a11*a11 + a22*a22) / 6.0
    + (a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]) / 3.0;
  /* solve x^3 - 3 p x  - 2 q = 0 */
  pr = sqrt(p);
  pr3 = p*pr;
  if (pr3 <= fabs(q)) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = (real) (m + pr);
      v[2] = (real) (m - 2.0 * pr);
    } else { /* phi = 0 */
      v[0] = (real) (m + 2.0 * pr);
      v[2] = v[1] = (real) (m - pr);
    }
  } else {
    double phi = acos(q/pr3)/3.0; /* 0 < phi < pi/3 */
    v[0] = (real) (m + 2.0 * pr * cos(phi));  /* largest */
    v[1] = (real) (m + 2.0 * pr * cos(phi - 2*M_PI/3)); /* second largest */
    v[2] = (real) (m + 2.0 * pr * cos(phi + 2*M_PI/3)); /* smallest */
#ifdef RV3_DEBUG
    {
      int i; double vi, y[3], dy[3];

      for (i = 0; i < 3; i++) { vi = v[i] - m; y[i] = vi*(vi*vi - 3*p) - 2 * q; dy[i] = 3*(vi*vi - p); }
      rv3_print(v,  "roots   ", "%26.14e", 1);
      dv3_print(y,  "residues", "%26.14e", 1);
      dv3_print(dy, "slope   ", "%26.14e", 1);
      printf("m %20.14f, q/pr3 %22.14e, pr %20.14f, phi %20.14f, %20.14f deg\n", m, q/pr3, pr, phi/M_PI*180, 120 - phi/M_PI*180);
    }
#endif
  }
  return v;
}



/* sort s to descending order, order u and v correspondingly */
INLINE void rv3_sort3(real s[3], real (*u)[3], real (*v)[3])
{
  real tmp;

  if (s[2] > s[1]) {
    tmp = s[1]; s[1] = s[2]; s[2] = tmp;
    if (u) rv3_swap(u[1], u[2]);
    if (v) rv3_swap(v[1], v[2]);
  }
  if (s[1] > s[0]) {
    tmp = s[0]; s[0] = s[1]; s[1] = tmp;
    if (u) rv3_swap(u[0], u[1]);
    if (v) rv3_swap(v[0], v[1]);
  }
  if (s[2] > s[1]) {
    tmp = s[1]; s[1] = s[2]; s[2] = tmp;
    if (u) rv3_swap(u[1], u[2]);
    if (v) rv3_swap(v[1], v[2]);
  }
}



/* return the pivot row for column c, row index starts from r0 */
INLINE int rm3_pivot_(real m[3][3], int r0, int c, real *max)
{
  int i, r = r0;
  real tmp;

  for (*max = fabs(m[r = r0][c]), i = r0 + 1; i < 3; i++)
    if ((tmp = fabs(m[i][c])) > *max) r = i, *max = tmp;
  return r;
}



/* solve matrix equation a x = 0, matrix 'a' is destroyed
 * solutions are saved as *row* vectors in 'x'
 * return the number of solutions */
INLINE int rm3_solvezero(real a[3][3], real (*x)[3], real tol)
{
  real max;
  int i, j, k, ns = 0;

  k = rm3_pivot_(a, 0, 0, &max); /* pivot for column 0 */
  if (max <= tol) { /* found the first eigenvector */
    rv3_make(x[ns++], 1, 0, 0);
    k = rm3_pivot_(a, 0, 1, &max); /* pivot for column 1 */
    if (max <= tol) {
      rv3_make(x[ns++], 0, 1, 0);
      k = rm3_pivot_(a, 0, 2, &max);
      if (max <= tol) rv3_make(x[ns++], 0, 0, 1);
    } else {
      if (k != 0) rv3_swap(a[0], a[k]);
      a[0][2] /= a[0][1]; /* normalize row 0, a[0][1] = 1; */
      for (i = 1; i < 3; i++) a[i][2] -= a[i][1]*a[0][2];
      rm3_pivot_(a, 1, 2, &max);
      if (max <= tol) rv3_makenorm(x[ns++], 0, a[0][2], -1);
    }
  } else {
    if (k != 0) rv3_swap(a[0], a[k]);
    a[0][1] /= a[0][0]; a[0][2] /= a[0][0]; /* normalize row 0, a[0][0] = 1 */
    for (i = 1; i < 3; i++)
      for (j = 1; j < 3; j++) /* a[i][0] = 0 now */
        a[i][j] -= a[i][0]*a[0][j];
    k = rm3_pivot_(a, 1, 1, &max); /* pivot for column 1 */
    if (max <= tol) { /* column 1 is empty */
      rv3_makenorm(x[ns++], -a[0][1], 1, 0);
      k = rm3_pivot_(a, 1, 2, &max);
      if (max <= tol) /* column 2 is empty too */
        rv3_makenorm(x[ns++], -a[0][2], 0, 1);
    } else {
      if (k != 1) rv3_swap(a[1], a[k]);
      a[1][2] /= a[1][1]; /* normalize row 1, a[1][1] = 1 */
      a[2][2] -= a[2][1]*a[1][2];
      if (fabs(a[2][2]) > tol) {
#ifdef RV3_DEBUG
        printf("a22 %g vs tol %g\n", a[2][2], tol);
#endif
        return 0; /* no solutions */
      }
      a[0][2] -= a[0][1]*a[1][2];
      rv3_makenorm(x[ns++], -a[0][2], -a[1][2], 1);
    }
  }
  return ns;
}


/* given an eigenvalue, return the corresponding eigenvectors */
INLINE int rm3_eigvecs(real (*vecs)[3], real mat[3][3], real val, real tol)
{
  real m[3][3];

  rm3_copy(m, mat); /* make a matrix */
  m[0][0] -= val; m[1][1] -= val; m[2][2] -= val;
  return rm3_solvezero(m, vecs, tol);
}



/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs'
 * ideally, eigenvalues should be sorted in magnitude-descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
INLINE rv3_t *rm3_eigsys(real v[3], real vecs[3][3], real mat[3][3], int nt)
{
  real vs[5][3], sq, tol, nv; /* for safety, vs needs 5 rows */
  int n = 0, nn, i;

  rm3_eigval(v, mat);
  for (sq = 0, i = 0; i < 3; i++) sq += rv3_sqr(mat[i]);
  /* errors of the eigenvalues from the cubic equation can reach sqrt(eps)
   * use a large tolerance */
  tol = 10.0 * sqrt(sq * DBL_EPSILON);

  for (nn = i = 0; i < 3; i++) {
    n = rm3_eigvecs(vs+nn, mat, v[nn], tol);
    if (n == 0) goto ERR;
    if ((nn += n) >= 3) break;
  }

  /* NOTE: make sure eigenvectors are orthogonal */
  rv3_normalize( rv3_cross(vs[2], vs[0], vs[1]) );
  rv3_normalize( rv3_cross(vs[1], vs[2], vs[0]) );

  rm3_copy(vecs, vs);
  for (i = 0; i < 3; i++) {
    nv = rv3_dot(rm3_mulvec(vs[i], mat, vecs[i]), vecs[i]);
    if (fabs(nv - v[i]) > tol) {
      fprintf(stderr, "corrupted eigenvalue i %d, %g vs. %g\n", i, nv, v[i]);
      goto ERR;
    }
#ifdef RV3_DEBUG
    printf("Eigenvalue: %22.14f vs %22.14f (corrected)\n", v[i], nv);
    rv3_print(vecs[i], "eigenvector i", "%20.12e", 1);
#endif
    v[i] = nv;
  }
#ifdef RV3_DEBUG
  printf("det(V) = %g\n", rm3_det(vecs));
#endif
  rv3_sort3(v, vecs, NULL);

  if (nt) return vecs; else return rm3_trans(vecs);
ERR:
  printf("fatal: bad eigenvalues, n %d, nn %d\n", n, nn);
  rm3_print(mat, "matrix", "%24.16e", 1);
  rv3_print(v, "eigenvalues", "%24.16e", 1);
  exit(1);
  return NULL;
}



/* SVD decomposition of a 3x3 matrix A = U S V^T */
INLINE void rm3_svd(real a[3][3], real u[3][3], real s[3], real v[3][3])
{
  int i, rank;
  real ata[3][3], us[3][3];

  /* 1. compute A^T A and its eigenvectors, which is V */
  rm3_tmul(ata, a, a);
  rm3_eigsys(s, v, ata, 1);
#ifdef RV3_DEBUG
  rv3_print(s, "S^2 ", "%22.14e", 1);
  rm3_print(ata, "A^T A ",  "%20.14f", 1);
  rm3_print(v, "V^T ",  "%20.14f", 1);
#endif

  /* 2. U^T = S^{-1} V^T A^T, and each row of U^T is an eigenvector
   * since eigenvectors are to be normalized, S^{-1} is unnecessary */
  if (s[0] <= 0.0) {
    rank = 0;
    rm3_copy(u, v);
  } else {
    double tol = 10. * sqrt(DBL_EPSILON);
    /* the test i = 1 + (s[1] > s[0]*tol) + (s[2] > s[0]*tol); */
    rm3_mult(u, v, a);
    for (i = 0; i < 3; i++) {
      rv3_copy(us[i], u[i]); /* save a copy of V^T A^T before normalizing it */
      s[i] = rv3_norm(u[i]);
      if (s[i] > 0.0) rv3_smul(u[i], 1.0/s[i]);
    }
    rank = 1;
    rank += (fabs(rv3_dot(u[0], u[1])) < tol && s[1] > tol);
    rank += (fabs(rv3_dot(u[0], u[2])) < tol && fabs(rv3_dot(u[1], u[2])) < tol && s[2] > tol);
#ifdef RV3_DEBUG
    rm3_print(u, "U^T ", "%22.14e", 1);
    rm3_print(us, "Us^T ", "%22.14e", 1);
    rv3_print(s, "S ", "%22.14e", 1);
    rm3_print(a, "A ",  "%20.14f", 1);
    printf("rank = %d, tol %g, det %g, u0.u0 %g, u0.u1 %g, u0.u2 %g, u1.u1 %g, u1.u2 %g, u2.u2 %g\n\n\n", rank, tol, rm3_det(u),
        rv3_sqr(u[0]), rv3_dot(u[0], u[1]), rv3_dot(u[0], u[2]), rv3_sqr(u[1]), rv3_dot(u[1], u[2]), rv3_sqr(u[2]));
#endif
    if (rank <= 2) {
      if (rank == 1) {
        real z[3] = {0, 0, 0}, w, tmp;

        w = fabs(u[0][i = 0]);
        if ((tmp = fabs(u[0][1])) < w) w = tmp, i = 1;
        if ((tmp = fabs(u[0][2])) < w) i = 2;
        z[i] = 1.0f; /* select the smallest element in u[0] as z */
        rv3_normalize( rv3_cross(u[1], z, u[0]) );
        s[1] = rv3_dot(u[1], us[1]); /* S = U^T (V^T A^T)^T is more accurate than sqrt(A^T A) */
        if (s[1] < 0) { s[1] = -s[1]; rv3_neg(u[1]); } /* make sure s[1] > 0 */
      }
      rv3_normalize( rv3_cross(u[2], u[0], u[1]) );
      s[2] = rv3_dot(u[2], us[2]);
      if (s[2] < 0) { s[2] = -s[2]; rv3_neg(u[2]); }
    }
    rv3_sort3(s, u, v);
#ifdef RV3_DEBUG
    printf("det(U) %g, det(V) %g\n", rm3_det(u), rm3_det(v));
#endif
  }
  rm3_trans(v);
  rm3_trans(u);
}



/* an old alias */
#define rotfit3 rv3_rmsd

/* least square fit from x to y after rotation/translation of the former
 * the best fit structure is saved to xf, if not NULL */
INLINE real rv3_rmsd(rv3_t *x, rv3_t *xf, rv3_t *y, const real *w, int n,
    real (*r)[3], real *t)
{
  int i;
  real wtot = 0, sq, dev = 0, dev0, detm;
  rv3_t xc, yc, xs, ys, sig, t_;
  real u[3][3], v[3][3], s[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, xy[3][3], r_[3][3];

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
    wtot = (real) n;
  } else {
    for (wtot = 0., i = 0; i < n; i++) {
      rv3_sinc(xc, x[i], w[i]);
      rv3_sinc(yc, y[i], w[i]);
      wtot += w[i];
    }
  }
  rv3_smul(xc, 1.f/wtot);
  rv3_smul(yc, 1.f/wtot);

  /* 2. compute 3x3 asymmetric covariance matrix S = (x-xc) (y-yc)^T */
  for (i = 0; i < n; i++) {
    rv3_diff(xs, x[i], xc); /* shift to the center avoid the translation */
    rv3_diff(ys, y[i], yc);
    rm3_vtv(xy, xs, ys);
    sq  = rv3_sqr(xs);
    sq += rv3_sqr(ys);
    if (w) {
      rm3_sinc(s, xy, w[i]);
      dev += w[i]*sq;
    } else {
      rm3_inc(s, xy);
      dev += sq; /* Tr(x^T x + y^T y) */
    }
  }
  dev0 = dev;

  /* 3. SVD decompose S = u sig v^T */
  rm3_svd(s, u, sig, v);

  /* 4. compute R = v u^T */
  rm3_mult(r, v, u);
  detm = rm3_det(r);

#define rmsd_dump_(title) { const char *rfmt = "%22.14e"; \
    printf("rmsd " title " fatal error: detm = %g, n = %d\n", detm, n); \
    rm3_print(r, "r", rfmt, 1); \
    printf("det(r) = %g\n", rm3_det(r)); \
    rm3_mult(r, u, v); rm3_print(r, "rx", rfmt, 1); \
    printf("det(rx) = %g\n", rm3_det(r)); \
    rm3_print(u, "u", rfmt, 1); \
    printf("det(u) = %g\n", rm3_det(u)); \
    rm3_print(v, "v", rfmt, 1); \
    printf("det(v) = %g\n", rm3_det(v)); \
    rm3_print(s, "s", rfmt, 1); \
    printf("det(s) = %g\n", rm3_det(s)); \
    rv3_print(sig, "sig", rfmt, 1); \
    exit(1); }
  if (fabs(fabs(detm) - 1) > 0.01) rmsd_dump_("bad svd");
  if (detm < 0) { /* to avoid a reflection */
    rm3_trans(u);
    rv3_neg(u[2]); /* flip the last eigenvector */
    rm3_mul(r, v, u);
    dev -= 2*(sig[0]+sig[1]-sig[2]);
    detm = rm3_det(r);
    if (fabs(fabs(detm) - 1) > 0.01) rmsd_dump_("bad inv.");
#undef rmsd_dump_
  } else {
    dev -= 2*(sig[0]+sig[1]+sig[2]); /* -2 Tr(R x y^T) */
  }
  if (dev < 0) dev = 0;
  rv3_diff(t, yc, rm3_mulvec(xs, r, xc)); /* t = yc - R xc */

  /* 5. compute the rotated structure */
  if (xf || dev < dev0*0.01) { /* if there's a large cancellation recompute the deviation */
    real xfit[3];
    for (dev = 0, i = 0; i < n; i++) {
      rv3_add(xfit, rm3_mulvec(xs, r, x[i]), t); /* xf = R x + t */
      sq = rv3_dist2(y[i], xfit);
      if (xf) rv3_copy(xf[i], xfit);
      dev +=  (w ? w[i]*sq : sq); /* recompute the deviation */
    }
  }
  return (real) sqrt(dev/wtot);
}



/* generate a random orthonormal (unitary) 3x3 matrix */
INLINE rv3_t *rm3_rnduni(real a[3][3])
{
  real dot;

  rv3_rnd(a[0], -.5f, 1.f);
  rv3_normalize(a[0]);

  rv3_rnd(a[1], -.5f, 1.f);
  /* normalize a[1] against a[0] */
  dot = rv3_dot(a[0], a[1]);
  rv3_sinc(a[1], a[0], -dot);
  rv3_normalize(a[1]);

  rv3_cross(a[2], a[0], a[1]);
  return a;
}



#endif

