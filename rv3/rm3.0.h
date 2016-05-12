#include "rv3.h"
#ifndef RM3_H__
#define RM3_H__
/* routines for 3x3 matrices */



#ifndef FM3_T
#define FM3_T fm3_t
typedef float fm3_t[3][3];
#endif

#ifndef DM3_T
#define DM3_T dm3_t
typedef double dm3_t[3][3];
#endif

#ifndef RM3_T
#define RM3_T rm3_t
typedef real rm3_t[3][3];
#endif


#define rm3_print(m, nm, fmt, nl) rm3_fprint(stdout, m, nm, fmt, nl)

INLINE void rm3_fprint(FILE *fp, real (*m)[3], const char *nm, const char *fmt, int nl)
{
  int i, j;

  if (nm) fprintf(fp, "%s:%c", nm, nl ? '\n' : ' ');
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      fprintf(fp, fmt, m[i][j]);
    fprintf(fp, "%s", nl ? "\n" : "; ");
  }
}



INLINE rv3_t *rm3_make(real rx[3][3], real a00, real a01, real a02,
    real a10, real a11, real a12, real a20, real a21, real a22)
{
  rv3_make(rx[0], a00, a01, a02);
  rv3_make(rx[1], a10, a11, a12);
  rv3_make(rx[2], a20, a21, a22);
  return rx;
}



#define rm3_makem(rx, x) rm3_make(rx, \
    (real) x[0][0], (real) x[0][1], (real) x[0][2], \
    (real) x[1][0], (real) x[1][1], (real) x[1][2], \
    (real) x[2][0], (real) x[2][1], (real) x[2][2])



/* zero matrix */
#define rm3_zero(x) rm3_makem(x, 0, 0, 0, 0, 0, 0, 0, 0, 0)



/* identity matrix */
#define rm3_one(x) rm3_makem(x, 1, 0, 0, 0, 1, 0, 0, 0, 1)



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
  rv3_smul2(a[0], v, u[0]);
  rv3_smul2(a[1], v, u[1]);
  rv3_smul2(a[2], v, u[2]);
  return a;
}



/* a += b */
INLINE rv3_t *rm3_inc(real a[3][3], real b[3][3])
{
  rv3_inc(a[0], b[0]);
  rv3_inc(a[1], b[1]);
  rv3_inc(a[2], b[2]);
  return a;
}



/* a += b*s */
INLINE rv3_t *rm3_sinc(real a[3][3], real b[3][3], real s)
{
  rv3_sinc(a[0], b[0], s);
  rv3_sinc(a[1], b[1], s);
  rv3_sinc(a[2], b[2], s);
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
INLINE real *rm3_mulvec(real * RESTRICT c, real a[3][3], const real *v)
{
  c[0] = rv3_dot(a[0], v);
  c[1] = rv3_dot(a[1], v);
  c[2] = rv3_dot(a[2], v);
  return c;
}



/* c = a^T v */
INLINE real *rm3_tmulvec(real * RESTRICT c, real a[3][3], const real *v)
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



/* rotate `x' around `u' by `ang', save the result to `y' */
INLINE real *rv3_rot(real * RESTRICT y, const real *x, const real *u, real ang)
{
  real m[3][3];

  rm3_mkrot(m, u, ang);
  rm3_mulvec(y, m, x);
  return y;
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
  pr3 = p * pr;
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
      int i;
      double vi, y[3], dy[3];

      for (i = 0; i < 3; i++) {
        vi = v[i] - m;
        y[i] = vi*(vi*vi - 3*p) - 2 * q;
        dy[i] = 3*(vi*vi - p);
      }
      rv3_print(v,  "roots   ", "%26.14e", 1);
      dv3_print(y,  "residues", "%26.14e", 1);
      dv3_print(dy, "slope   ", "%26.14e", 1);
      printf("m %20.14f, q/pr3 %22.14e, pr %20.14f, phi %20.14f, %20.14f deg\n", m, q/pr3, pr, phi/M_PI*180, 120 - phi/M_PI*180);
    }
#endif
  }
  return v;
}



/* sort `s' to descending order, order `u' and `v' correspondingly */
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



/* return the pivot row r and column c, starting from (r0, c0)
 * cmap[r0] registers the actual column index */
INLINE double rm3_pivot_(real m[3][3], int r0, int cmap[])
{
  int i, j, r, c;
  double tmp, max;
  real t;

  /* 1. find the pivot row and column */
  max = -1;
  for (j = r0; j < 3; j++)
    for (i = r0; i < 3; i++)
      if ((tmp = fabs(m[i][j])) > max)
        r = i, c = j, max = tmp;
#ifdef RV3_DEBUG
  printf("r0 %d, r %d, c %d, max %g\n", r0, r, c, max);
#endif

  /* 2. put the pivot to the left-top corner */
  /* swap rows r and r0, which doesn't affect the solution */
  if (r != r0) rv3_swap(m[r], m[r0]);

  if (c != r0) { /* swap columns c and r0 */
    for (i = 0; i < 3; i++) /* must be from row 0 */
      t = m[i][c], m[i][c] = m[i][r0], m[i][r0] = t;
    i = cmap[c], cmap[c] = cmap[r0], cmap[r0] = i;
  }
  return max;
}



/* Solve matrix equation a x = 0 by Gaussian elimination (full-pivot)
 * The matrix 'a' is destroyed, solutions are saved as *row* vectors in 'x'
 * return the number of solutions
 * Note that this routine assumes at least one solution */
INLINE int rm3_solvezero(real a[3][3], real (*x)[3], real tol)
{
  double max;
  int cmap[3] = {0, 1, 2};
  real a00;

#ifdef RV3_DEBUG
  rm3_print(a, "input matrix", "%20.12f", 1);
#endif
  max = rm3_pivot_(a, 0, cmap); /* pivot for column 0 */
#ifdef RV3_DEBUG
  rm3_print(a, "after the first pivot", "%20.12f", 1);
  printf("pivot %d, max %f\n", cmap[0], max);
#endif

  if ( max <= tol ) { /* matrix is zero */
    rv3_make(x[0], 1, 0, 0);
    rv3_make(x[1], 0, 1, 0);
    rv3_make(x[2], 0, 0, 1);
    return 3;
  }

  /* normalize row 0 such that a[0][0] = 1 */
  a00 = a[0][0];
  a[0][1] /= a00;
  a[0][2] /= a00;
#ifdef RV3_DEBUG
  rm3_print(a, "normalized a", "%20.12f", 1);
#endif
  /* gaussian elimination */
  a[1][1] -= a[1][0] * a[0][1];
  a[1][2] -= a[1][0] * a[0][2];
  /* we don't really need a[2], just for pivoting */
  a[2][1] -= a[2][0] * a[0][1];
  a[2][2] -= a[2][0] * a[0][2];
#ifdef RV3_DEBUG
  a[0][0] = 1;
  a[1][0] = 0;
  a[2][0] = 0;
  rm3_print(a, "after first gaussian elimination", "%20.12f", 1);
#endif

  max = rm3_pivot_(a, 1, cmap); /* pivot for column 1 */
#ifdef RV3_DEBUG
  rm3_print(a, "after second pivot", "%20.12f", 1);
  printf("second pivot %d, max %g, tol %g\n",
      cmap[1], max, tol);
#endif

  if ( max <= tol ) { /* zero */
    a[1][ cmap[0] ] = 1; /* a[0][0] */
    a[1][ cmap[1] ] = a[0][1];
    a[1][ cmap[2] ] = a[0][2];
    /* the first eigenvector x[0] is normal to a[1]
     * since a[0][0] = 1 is the largest component, the
     * following vector is not vanishing */
    x[0][ cmap[0] ] = -a[0][1];
    x[0][ cmap[1] ] = 1; /* a[0][0] */
    x[0][ cmap[2] ] = 0;
    rv3_normalize( x[0] );
    /* the second eigenvector x[1] is normal to a[1] and x[0] */
    rv3_normalize( rv3_cross(x[1], a[1], x[0]) );
#ifdef RV3_DEBUG
    rv3_print(a[1], "degen2. vector 0", "%20.12f", 1);
    rv3_print(x[0], "degen2. vector 1", "%20.12f", 1);
    rv3_print(x[1], "degen2. vector 2", "%20.12f", 1);
#endif
    return 2;
  }

  x[0][ cmap[1] ] = a[1][2];
  x[0][ cmap[2] ] = -a[1][1];
  /* this component is used to make x[0] normalize to a[0] */
  x[0][ cmap[0] ] = - a[0][1] * a[1][2] - a[0][2] * (-a[1][1]);
  rv3_normalize(x[0]);
#ifdef RV3_DEBUG
  printf("cmap %d, %d, %d\n", cmap[0], cmap[1], cmap[2]);
  rv3_print(x[0], "eigenvector", "%20.12f", 1);
#endif
  return 1;
}



float fm3_eigtol = 10 * FLT_EPSILON;
double dm3_eigtol = 10 * DBL_EPSILON;
real rm3_eigtol = 10 * DBL_EPSILON;



#define rm3_eigvecs(vecs, mat, val) \
  rm3_eigvecs_(vecs, mat, val, rm3_eigtol)



/* compute the eigenvector of a given eigvalue `val`
 * no matter how small `reltol` is, this routine
 * returns at least one eigenvector, whereas more
 * eigenvectors are possible with a larger `reltol` */
INLINE int rm3_eigvecs_(real (*vecs)[3], real mat[3][3], real val,
    real reltol)
{
  real m[3][3], max = 0, x, tol;
  int i, j;

  /* find the maximal element of m */
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ )
      if ( (x = (real) fabs(mat[i][j])) > max )
        max = x;
  tol = max * reltol;

  rm3_copy(m, mat); /* make a matrix */
  m[0][0] -= val;
  m[1][1] -= val;
  m[2][2] -= val;
#ifdef RV3_DEBUG
  printf("\n\nsolving eigenvector(s) for the eigenvalue %g\n", val);
#endif
  return rm3_solvezero(m, vecs, tol);
}



#define rm3_eigsys(vals, vecs, mat, nt) \
  rm3_eigsys_(vals, vecs, mat, nt, rm3_eigtol)

/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs'
 * ideally, eigenvalues should be sorted in magnitude-descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
INLINE rv3_t *rm3_eigsys_(real vals[3], real vecs[3][3], real mat[3][3],
    int nt, real reltol)
{
  real vs[3][3] = {{0}}, x;
  int n;

  /* 1. compute the eigenvalues from solving the cubic equation
   * the eigenvectors are descending */
  rm3_eigval(vals, mat);

  /* 2. solve the eigenvector of vals[0]
   * this routine returns at least one solution */
  n = rm3_eigvecs_(vecs, mat, vals[0], reltol);

  if ( n == 1 ) { /* we only got one vector */
    /* 3. solve the eigenvector of vals[1] */
    /* swap the last two eigenvalues to avoid degeneracy
     * with vals[0] and vals[1] well separated, degeneracy is less likely */
    x = vals[1], vals[1] = vals[2], vals[2] = x;
    n = rm3_eigvecs_(vs, mat, vals[1], reltol);

    /* check if we get a degeneracy, i.e., the eigenvector is
     * essentially the same as the previous one */
    if ( (x = (real) fabs(rv3_dot(vs[0], vecs[0]))) > 0.5 ) {
      /* This happens only if the largest eigenvalue vals[0]
       * and smallest vals[1] are essentially the same,
       * and the tolerance is set too small
       * which means all eigenvectors are degenerate
       * 0.5 is a very loose threshold */
      rv3_make(vecs[0], 1, 0, 0);
      rv3_make(vecs[1], 0, 1, 0);
      rv3_make(vecs[2], 0, 0, 1);
    } else {
      /* compute the last eigenvector */
      rv3_copy( vecs[1], vs[0] );
      rv3_normalize( rv3_cross(vecs[2], vecs[0], vecs[1]) );
    }
  } else if ( n == 2 ) {
    /* we already have two vectors, deduce the last */
    rv3_normalize( rv3_cross(vecs[2], vecs[0], vecs[1]) );
  } /* if n == 3, we are done */

#ifdef RV3_DEBUG
  {
    int i;
    for (i = 0; i < 3; i++) {
      double tol = 1e-6, sq, nv;
      int j;
      nv = rv3_dot(rm3_mulvec(vs[i], mat, vecs[i]), vecs[i]);
      if (sizeof(real) == sizeof(float)) tol = 1e-4;
      for (sq = 0, j = 0; j < 3; j++) sq += rv3_sqr(mat[j]);
      tol *= sqrt(sq)/3;
      if (fabs(nv - vals[i]) > tol) {
        fprintf(stderr, "corrupted eigenvalue i %d, %g vs. %g (%g, %g, %g)\n",
            i, nv, vals[i], vals[0], vals[1], vals[2]);
        rm3_print(mat, "matrix", "%20.14f", 1);
        rv3_print(vals, "eigenvalues  ", "%20.14f", 1);
        rv3_print(vecs[0], "eigenvector 0", "%20.14f", 1);
        rv3_print(vecs[1], "eigenvector 1", "%20.14f", 1);
        rv3_print(vecs[2], "eigenvector 2", "%20.14f", 1);
        return NULL;
      }
      vals[i] = nv;
    }
    printf("det(V) = %g\n", rm3_det(vecs));
  }
#endif
  rv3_sort3(vals, vecs, NULL);

  return nt ? vecs : rm3_trans(vecs);
}



/* SVD decomposition of a 3x3 matrix A = U S V^T */
INLINE void rm3_svd(real a[3][3], real u[3][3], real s[3], real v[3][3])
{
  int i, rank;
  real ata[3][3], us[3][3];

  /* A^T A = V S^2 V^T, so (A^T A) V = V S^2 */

  /* 1. compute A^T A and its eigenvectors, which is V */
  rm3_tmul(ata, a, a);
  rm3_eigsys(s, v, ata, 1);
#ifdef RV3_DEBUG
  rv3_print(s,   "S^2 ",   "%22.12e", 1);
  rm3_print(ata, "A^T A ", "%22.14f", 1);
  rm3_print(v,   "V^T ",   "%22.14f", 1);
  printf("det(V) = %g\n\n\n", rm3_det(v));
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
      if (s[i] > 0) rv3_smul(u[i], 1/s[i]);
    }
    rank = 1;
    rank += (fabs(rv3_dot(u[0], u[1])) < tol && s[1] > tol);
    rank += (fabs(rv3_dot(u[0], u[2])) < tol && fabs(rv3_dot(u[1], u[2])) < tol && s[2] > tol);
#ifdef RV3_DEBUG
    rm3_print(u,  "U^T ",  "%22.15f", 1);
    rm3_print(us, "Us^T ", "%22.15f", 1);
    rv3_print(s,  "S ",    "%22.15f", 1);
    rm3_print(a,  "A ",    "%22.15f", 1);
    printf("rank = %d, tol %g, det %g, u0.u0 %g, u0.u1 %g, u0.u2 %g, u1.u1 %g, u1.u2 %g, u2.u2 %g\n\n\n", rank, tol, rm3_det(u),
        rv3_sqr(u[0]), rv3_dot(u[0], u[1]), rv3_dot(u[0], u[2]), rv3_sqr(u[1]), rv3_dot(u[1], u[2]), rv3_sqr(u[2]));
#endif
    if (rank <= 2) {
      if (rank == 1) {
        real z[3] = {0, 0, 0}, w, tmp;

        w = (real) fabs(u[0][i = 0]);
        if ((tmp = (real) fabs(u[0][1])) < w) w = tmp, i = 1;
        if ((tmp = (real) fabs(u[0][2])) < w) i = 2;
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

/* Fit x to y by rotation and translation of the `x'
 * If `refl', reflection can also be used.
 * The best-fit structure is saved to `xf', if not NULL */
INLINE real rv3_rmsd(rv3_t * RESTRICT x, rv3_t * RESTRICT xf,
    rv3_t * RESTRICT y, const real *w, int n, int refl,
    real (* RESTRICT r)[3], real * RESTRICT t)
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
  rv3_diff(t, yc, xc);

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

#define rmsd_dump_(title) { const char *rfmt = "%22.14f"; \
    printf("rmsd [" title "], fatal error: detm = %g, n = %d\n", detm, n); \
    rm3_print(s, "s = (x - xc) (y - yc)^T", rfmt, 1); \
    printf("det(s) = %g\n", rm3_det(s)); \
    rm3_print(u, "u", rfmt, 1); \
    printf("det(u) = %g\n", rm3_det(u)); \
    rm3_print(v, "v", rfmt, 1); \
    printf("det(v) = %g\n", rm3_det(v)); \
    rv3_print(sig, "sig", rfmt, 1); \
    rm3_print(r, "r = v.u (rotation matrix)", rfmt, 1); \
    printf("det(r) = %g\n", rm3_det(r)); \
    rm3_mult(r, u, v); rm3_print(r, "r' = u.v", rfmt, 1); \
    printf("det(r') = %g\n", rm3_det(r)); \
    exit(1); }
  if (fabs(fabs(detm) - 1) > 0.01) {
    fprintf(stderr, "detm: %g\n", detm);
    rmsd_dump_("bad svd");
  }
  if (detm < 0 && !refl) { /* to avoid a reflection */
    rm3_trans(u);
    rv3_neg(u[2]); /* flip the last eigenvector */
    rm3_mul(r, v, u);
    dev -= 2*(sig[0] + sig[1] - sig[2]);
    detm = rm3_det(r);
    if (fabs(fabs(detm) - 1) > 0.01) rmsd_dump_("bad inv.");
#undef rmsd_dump_
  } else {
    dev -= 2 * (sig[0] + sig[1] + sig[2]); /* -2 Tr(R x y^T) */
  }
  if (dev < 0) dev = 0;

  /* 5. compute the rotated structure */
  if (xf || dev < dev0*0.01) { /* if there's a large cancellation recompute the deviation */
    real xfi[3];

    for (dev = 0, i = 0; i < n; i++) {
      rv3_diff(xs, x[i], xc);
      rm3_mulvec(ys, r, xs);
      rv3_add(xfi, ys, yc); /* xfi = R (x - xc) + yc */
      sq = rv3_dist2(y[i], xfi);
      if (xf) rv3_copy(xf[i], xfi);
      dev +=  (w ? w[i]*sq : sq); /* recompute the deviation */
    }
  }
  return (real) sqrt(dev/wtot);
}



/* generate a random orthonormal (unitary) 3x3 matrix */
INLINE rv3_t *rm3_randuni(real a[3][3])
{
  real dot;

  rv3_randdir0(a[0]);

  rv3_randunif(a[1], -1, 1);
  /* component of a[1] normal to a[0] */
  dot = rv3_dot(a[0], a[1]);
  rv3_normalize( rv3_sinc(a[1], a[0], -dot) );

  rv3_cross(a[2], a[0], a[1]);
  return a;
}



#endif

