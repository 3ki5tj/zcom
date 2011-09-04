#ifndef RV2_H__
#define RV2_H__

#ifdef HAVE_REAL
  #ifndef ZCHAVEREAL
  #define ZCHAVEREAL HAVE_REAL
  #endif
#endif

#ifndef ZCHAVEREAL
  #define ZCHAVEREAL 1
  typedef double real;
#endif

#ifndef RV2_T
#define RV2_T rv2_t
  typedef real rv2_t[2];
  typedef const real crv2_t[2];
#endif

#include <math.h>
#include <string.h>

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */
ZCINLINE real *rv2_make(real *x, real a, real b) { x[0] = a; x[1] = b; return x; }
ZCINLINE real *rv2_zero(real *x) { return rv2_make(x, 0, 0); }
ZCINLINE real *rv2_copy(real *x, const real *src) { x[0] = src[0]; x[1] = src[1]; return x; }
/* use macro to avoid const qualifier of src */
#define rv2_ncopy(x, src, n) memcpy(x, src, n*sizeof(x[0]))

ZCINLINE real rv2_sqr(const real *x) { return x[0]*x[0]+x[1]*x[1]; }
ZCINLINE real rv2_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]); }

ZCINLINE real rv2_dot(const real *x, const real *y) { return x[0]*y[0]+x[1]*y[1]; }

ZCINLINE real rv2_cross(const real *x, const real *y)
{
  return x[0]*y[1]-x[1]*y[0];
}

ZCINLINE real *rv2_neg(real *x)
{
  x[0] -= x[0];
  x[1] -= x[1];
  return x;
}

ZCINLINE real *rv2_neg2(real *nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  return nx;
}

ZCINLINE real *rv2_inc(real *x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  return x;
}

ZCINLINE real *rv2_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  return x;
}

ZCINLINE real *rv2_sinc(real *x, const real *dx, real s)
{
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  return x;
}

ZCINLINE real *rv2_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  return x;
}

ZCINLINE real *rv2_smul2(real *y, const real *x, real s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  return y;
}

ZCINLINE real *rv2_normalize(real *x)
{
  real r = rv2_norm(x);
  if (r > 0.f) rv2_smul(x, 1./r);
  return x;
}

/* for in-place difference use rv2_dec */
ZCINLINE real *rv2_diff(real *diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  return diff;
}

/* distance^2 between a and b */
ZCINLINE real rv2_dist2(const real *a, const real *b) 
{
  real d[2]; 
  return rv2_sqr(rv2_diff(d, a, b));
}

/* distance between a and b */
ZCINLINE real rv2_dist(const real *a, const real *b) 
{
  return (real) sqrt(rv2_dist2(a, b));
}

/* sum = a+b, for in-place addition use rv2_inc */
ZCINLINE real *rv2_add(real *sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv2_nadd(real *sum, const real *a, const real *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  return sum;
}

ZCINLINE real *rv2_lincomb2(real *sum, const real *a, const real *b, real s1, real s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  return sum;
}

/* consine of the angle of x1-x2-x3 */
ZCINLINE real rv2_cosang(const real *x1, const real *x2, const real *x3, 
    real *g1, real *g2, real *g3)
{
  real a[2], b[2], ra, rb, dot;

  ra = rv2_norm(rv2_diff(a, x1, x2));
  rv2_smul(a, 1.f/ra);
  rb = rv2_norm(rv2_diff(b, x3, x2));
  rv2_smul(b, 1.f/rb);
  dot = rv2_dot(a, b);
  if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
  if (g1) {
    rv2_lincomb2(g1, b, a, 1.f/ra, -dot/ra);
    rv2_lincomb2(g3, a, b, 1.f/rb, -dot/rb);
    rv2_nadd(g2, g1, g3);
  }
  return dot;
}

/* angle and gradients of x1-x2-x3 */
ZCINLINE real rv2_ang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real dot, sn;
  dot = rv2_cosang(x1, x2, x3, g1, g2, g3);
  sn = (real) sqrt(1 - dot*dot);
  if (sn > 1e-7) sn = -1/sn; else sn = 0.;
  if (g1) {
    rv2_smul(g1, sn);
    rv2_smul(g2, sn);
    rv2_smul(g3, sn);
  }
  return (real) acos(dot);
}

/* vertical distance from x to line a-b */
ZCINLINE real rv2_vdist(const real *x, const real *a, const real *b)
{
  real nm[2], d[2], dot;

  rv2_diff(d, x, a);
  rv2_normalize(rv2_diff(nm, a, b));
  dot = rv2_dot(d, nm);
  return rv2_norm(rv2_sinc(d, nm, -dot));
}

/* determinant of a 2x2 matrix */
ZCINLINE real mat2_det(real a[2][2])
{
  return a[0][0]*a[1][1] - a[0][1]*a[1][0];
}

/* inverse matrix b = a^(-1) */
ZCINLINE void mat2_inv(real b[2][2], real a[2][2])
{
  real det = mat2_det(a);
  if (fabs(det) < 1e-30) det = (det < 0) ? -1e-30f: 1e-30f;
  b[0][0] =  a[1][1]/det;
  b[0][1] = -a[0][1]/det;
  b[1][0] = -a[1][0]/det;
  b[1][1] =  a[0][0]/det;
}

#endif /* RV2_H__ */

