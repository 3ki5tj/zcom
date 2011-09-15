#define ZCINLINE __inline static
#define ZCRESTRICT __restrict 
#include "def.h"
#ifndef RV3_H__
#define RV3_H__

#ifndef RV3_T
#define RV3_T rv3_t
  typedef real rv3_t[3];
  typedef const real crv3_t[3];
  typedef real mat3_t[3][3];
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */

ZCINLINE real *rv3_make(real *x, real a, real b, real c)
  { x[0] = a; x[1] = b; x[2] = c; return x; }
ZCINLINE real *rv3_zero(real *x) { return rv3_make(x, 0, 0, 0); }
ZCINLINE real *rv3_copy(real *x, const real *src)
  { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; return x; }
/* use macro to avoid const qualifier of src */
#define rv3_ncopy(x, src, n) memcpy(x, src, n*sizeof(x[0]))

ZCINLINE real rv3_sqr (const real *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
ZCINLINE real rv3_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

/* if x == y, try to use sqr */
ZCINLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

ZCINLINE real *rv3_cross(real *ZCRESTRICT z, const real *x, const real *y)
{
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
  return z;
}

ZCINLINE real *rv3_neg(real *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  x[2] = -x[2];
  return x;
}

ZCINLINE real *rv3_neg2(real *nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  nx[2] = -x[2];
  return nx;
}

ZCINLINE real *rv3_inc(real *x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  x[2] += dx[2];
  return x;
}

ZCINLINE real *rv3_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  x[2] -= dx[2];
  return x;
}

ZCINLINE real *rv3_sinc(real *x, const real *dx, real s)
{
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  x[2] += s*dx[2];
  return x;
}

ZCINLINE real *rv3_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}

/* if y == x, just use smul */
ZCINLINE real *rv3_smul2(real * ZCRESTRICT y, const real *x, real s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  y[2] = x[2]*s;
  return y;
}

ZCINLINE real *rv3_normalize(real *x)
{
  real r = rv3_norm(x);
  if (r > 0.0) rv3_smul(x, 1.f/r);
  return x;
}

/* for in-place difference use rv3_dec */
ZCINLINE real *rv3_diff(real * ZCRESTRICT diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  diff[2] = a[2]-b[2];
  return diff;
}

/* distance^2 between a and b */
ZCINLINE real rv3_dist2(const real *a, const real *b) 
{
  real d[3]; 
  return rv3_sqr(rv3_diff(d, a, b));
}

/* distance between a and b */
ZCINLINE real rv3_dist(const real *a, const real *b) 
{
  return (real) sqrt(rv3_dist2(a, b));
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv3_add(real * ZCRESTRICT sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  sum[2] = a[2]+b[2];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv3_nadd(real *sum, const real *a, const real *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  sum[2] = -a[2]-b[2];
  return sum;
}

ZCINLINE real *rv3_lincomb2(real *sum, const real *a, const real *b, real s1, real s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  sum[2] = a[2]*s1+b[2]*s2;
  return sum;
}

/* angle and gradients of cos(x1-x2-x3) */
ZCINLINE real rv3_cosang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real a[3], b[3], ra, rb, dot;

  ra = rv3_norm(rv3_diff(a, x1, x2));
  rv3_smul(a, 1.f/ra);
  rb = rv3_norm(rv3_diff(b, x3, x2));
  rv3_smul(b, 1.f/rb);
  dot = rv3_dot(a, b);
  if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
  if (g1) {
    rv3_lincomb2(g1, b, a, 1.f/ra, -dot/ra);
    rv3_lincomb2(g3, a, b, 1.f/rb, -dot/rb);
    rv3_nadd(g2, g1, g3);
  }
  return dot;
}

/* angle and gradients of x1-x2-x3 */
ZCINLINE real rv3_ang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real dot, sn;
  
  dot = rv3_cosang(x1, x2, x3, g1, g2, g3);
  sn = (real) sqrt(1 - dot*dot);
  if (sn < 1e-7) sn = 1; else sn = -1.f/sn;
  if (g1) {
    rv3_smul(g1, sn);
    rv3_smul(g2, sn);
    rv3_smul(g3, sn);
  }
  return (real) acos(dot);
}

/* vertical distance from x to line a-b */
ZCINLINE real rv3_vdist(const real *x, const real *a, const real *b)
{
  real nm[3], d[3], dot;

  rv3_diff(d, x, a);
  rv3_normalize(rv3_diff(nm, a, b));
  dot = rv3_dot(d, nm);
  return rv3_norm(rv3_sinc(d, nm, -dot));
}

/* signed distance from x to the plane extended by a, b, c */
ZCINLINE real rv3_vpdist(const real *x, const real *a, const real *b, const real *c)
{
  real u[3], v[3], m[3];

  rv3_diff(u, b, a);
  rv3_diff(v, c, b);
  rv3_normalize(rv3_cross(m, u, v));
  rv3_diff(u, x, a);
  return rv3_dot(u, m);
}

/* transpose */
ZCINLINE rv3_t *mat3_trans(real a[3][3]) 
{
  real x;
  x = a[0][1], a[0][1] = a[1][0], a[1][0] = x;
  x = a[0][2], a[0][2] = a[2][0], a[2][0] = x;
  x = a[2][1], a[2][1] = a[1][2], a[1][2] = x;
  return a;
}

/* a = u^T v */
ZCINLINE rv3_t *mat3_vtv(real a[3][3], const real *u, const real *v)
{
  a[0][0] = u[0]*v[0];
  a[0][1] = u[0]*v[1];
  a[0][2] = u[0]*v[2];
  a[1][0] = u[1]*v[0];
  a[1][1] = u[1]*v[1];
  a[1][2] = u[1]*v[2];
  a[2][0] = u[2]*v[0];
  a[2][1] = u[2]*v[1];
  a[2][2] = u[2]*v[2];
  return a;
}

/* a += b */
ZCINLINE rv3_t *mat3_inc(real a[3][3], real b[3][3])
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
ZCINLINE rv3_t *mat3_sinc(real a[3][3], real b[3][3], real s)
{
  a[0][0] += b[0][0]*s;
  a[0][1] += b[0][1]*s;
  a[0][2] += b[0][2]*s;
  a[1][0] += b[1][0]*s;
  a[1][1] += b[1][1]*s;
  a[1][2] += b[1][2]*s;
  a[2][0] += b[2][0]*s;
  a[2][1] += b[2][1]*s;
  a[2][2] += b[2][2]*s;
  return a;
}

/* c = a b */
ZCINLINE rv3_t *mat3_mul(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
  return c;
}

/* c = a b^T */
ZCINLINE rv3_t *mat3_mult(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = rv3_dot(a[i], b[j]);
  return c;
}

/* c = a v */
#define rv3_rot(v1, rot, v) mat3_mulvec(v1, rot, v)
ZCINLINE real *mat3_mulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0]+a[0][1]*v[1]+a[0][2]*v[2];
  c[1] = a[1][0]*v[0]+a[1][1]*v[1]+a[1][2]*v[2];
  c[2] = a[2][0]*v[0]+a[2][1]*v[1]+a[2][2]*v[2];
  return c;
}

/* c = a^T v */
ZCINLINE real *mat3_multvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0]+a[1][0]*v[1]+a[2][0]*v[2];
  c[1] = a[0][1]*v[0]+a[1][1]*v[1]+a[2][1]*v[2];
  c[2] = a[0][2]*v[0]+a[1][2]*v[1]+a[2][2]*v[2];
  return c;
}

/* determinant of a 3x3 matrix */
ZCINLINE real mat3_det(real a[3][3])
{
  return a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[0][2]*a[1][2]+a[1][0]*a[2][0]*a[2][1]
    - (a[0][0]*a[1][2]*a[2][1]+a[1][1]*a[0][2]*a[2][0]+a[2][2]*a[0][1]*a[1][0]);
}

/* inverse matrix b = a^(-1) */
ZCINLINE rv3_t *mat3_inv(real b[3][3], real a[3][3])
{
  real dt = mat3_det(a);
  if (fabs(dt) < 1e-30) dt = (dt < 0) ? -1e-30f: 1e-30f;
  b[0][0] = (a[1][1]*a[2][2] - a[1][2]*a[2][1])/dt;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/dt;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/dt;
  b[1][0] = (a[1][2]*a[2][0] - a[1][0]*a[2][2])/dt;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/dt;
  b[2][0] = (a[1][0]*a[2][1] - a[2][0]*a[1][1])/dt;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/dt;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/dt;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/dt;
  return b;
}

/* eigenvalues of a 3x3 matrix */
ZCINLINE real *mat3_eigval(real v[3], real a[3][3])
{
  real m, p, q, cphi, sphi, pr, pr3;

  m = (real)((a[0][0]+a[1][1]+a[2][2])/3.);
  a[0][0] -= m;
  a[1][1] -= m;
  a[2][2] -= m;
  q = .5f*mat3_det(a);
  p = ((a[0][0]*a[0][0]+a[1][1]*a[1][1]+a[2][2]*a[2][2]) +
    2.f*(a[0][1]*a[1][0]+a[2][0]*a[0][2]+a[1][2]*a[2][1]))/6.f;
  pr = (real)sqrt(p);
  pr3 = p*pr;
  if (pr3 <= fabs(q)) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = m + pr;
      v[2] = m - 2.f*pr;
    } else { /* phi = 0 */
      v[0] = m + 2.f*pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    double phi = acos(q/pr3)/3.f; /* 0 < phi < pi/3 */
    cphi = (real)cos(phi);
    sphi = (real)(sin(phi)*1.7320508075688772);
    v[0] = m + 2.f*pr*cphi;  /* cos(phi), largest */
    v[1] = m - pr*(cphi-sphi); /* cos(phi-2*pi/3), second largest */
    v[2] = m - pr*(cphi+sphi); /* cos(phi+2*pi/3), smallest */
  }
  a[0][0] += m;
  a[1][1] += m;
  a[2][2] += m;  
  return v;
}

/* given matrix a and eigenvalue lm, return eigenvector */
ZCINLINE real *mat3_eigvec(real vec[3], real m[3][3], real val)
{
  double a = m[0][0]-val, b = m[1][1]-val, c, d = m[0][1], e = m[0][2], f = m[1][2];
  double dt, tol = 1e-12;

  vec[2] = 1.f;
  if (fabs(dt = a*b - d*d) > tol) { /* use row 0 and 1 */
    rv3_make(vec, (real)((d*f-b*e)/dt), (real)((e*d-a*f)/dt), 1);
    return rv3_normalize(vec);
  }
  c = m[2][2] - val;
  if (fabs(dt = a*f - e*d) > tol) { /* row 1 and 2 */
    rv3_make(vec, (real)((d*c-e*f)/dt), (real)((e*e-a*c)/dt), 1);
    return rv3_normalize(vec);
  }
  if ((dt = sqrt(a*a+d*d)) > tol) { /* three-row-degenerate */
    rv3_make(vec, (real)(d/dt), (real)(-a/dt), 0);
  } else {
    rv3_make(vec, 1, 0, 0);
  }
  return vec;
}

/* compute eigenvectors for the eigenvalues */
ZCINLINE rv3_t *mat3_eigvecs(real vecs[3][3], real mat[3][3], real v[3], int t)
{
  const double tol = 1e-12;
  double v0 = fabs(v[0]), v1 = fabs(v[1]), v2 = fabs(v[2]);
  
  mat3_eigvec(vecs[0], mat, v[0]);
  if ( fabs(v[0] - v[1]) > tol*(v0 + v1) ) {
    mat3_eigvec(vecs[1], mat, v[1]);
    rv3_cross(vecs[2], vecs[0], vecs[1]);
  } else if ( fabs(v[2] - v[1]) > tol*(v1 + v2) ) {
    mat3_eigvec(vecs[2], mat, v[2]);
    rv3_cross(vecs[1], vecs[2], vecs[0]);
  } else {
    rv3_make(vecs[1], 0, 1, 0);
    rv3_make(vecs[2], 0, 0, 1);
  }
  /* transpose the matrix */
  if (t) return vecs;
  else return mat3_trans(vecs);
}

/* SVD decomposition of a 3x3 matrix a = u s v */
ZCINLINE int mat3_svd(real a[3][3], real u[3][3], real s[3], real v[3][3])
{
  int i, j;
  real ata[3][3], z[3];
  const double tol = 1e-12;

  /* 1. compute A^T A and its eigenvectors */
  for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) {
      ata[i][j] = a[0][i]*a[0][j] + a[1][i]*a[1][j] + a[2][i]*a[2][j];
      if (i != j) ata[j][i] = ata[i][j];
    }
  mat3_eigval(s, ata);
  for (i = 0; i < 3; i++) if (s[i] < 0) s[i] = 0;
  mat3_eigvecs(v, ata, s, 1); /* get V^T */

  /* 2. U = A V S^-1, or U^T = S^{-1}T V^T A^T */
  j = (s[0] > tol) + (s[1] > tol) + (s[2] > tol);
  if (j >= 2) {
    mat3_mult(u, v, a);
    if (j == 2) rv3_cross(u[2], u[0], u[1]); /* fix the last */
  } else if (j == 1) {
    mat3_multvec(u[0], a, v[0]);
    rv3_zero(z);
    /* choose z[i] such that z X u[0] != 0 */
    i = (u[0][0]*u[0][0] < u[0][1]*u[0][1]) ? 0 : 1;
    rv3_cross(u[1], z, u[0]);
    rv3_cross(u[2], u[0], u[1]);
  } else { /* use u */
    for (i = 0; i < 3; i++) rv3_copy(u[i], v[i]);
  }
  for (i = 0; i < 3; i++) rv3_normalize(u[i]);
  for (i = 0; i < 3; i++) s[i] = (real)sqrt(s[i]);
  mat3_trans(u);
  mat3_trans(v);
  return 0;
} 

#define rv3_print(r, nm, fmt, nl) rv3_fprint(stdout, r, nm, fmt, nl)
ZCINLINE void rv3_fprint(FILE *fp, real r[3], const char *nm,
    const char *fmt, int nl)
{
  int i;
  if (nm) fprintf(fp, "%s: ", nm);
  for (i = 0; i < 3; i++)
    fprintf(fp, fmt, r[i], nl);
  fprintf(fp, "%c", (nl ? '\n' : ';'));
}

#define mat3_print(r, nm, fmt, nl) mat3_fprint(stdout, r, nm, fmt, nl)
ZCINLINE void mat3_fprint(FILE *fp, real r[3][3], const char *nm,
    const char *fmt, int nl)
{
  int i, j;
  if (nm) fprintf(fp, "%s:%c", nm, (nl ? '\n' : ' '));
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      fprintf(fp, fmt, r[i][j], nl);
    }
    fprintf(fp, "%s", (nl ? "\n" : "; "));
  }
}

#endif /* RV3_H__ */

