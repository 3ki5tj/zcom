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

#include <math.h>

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */

ZCINLINE void rv2_zero(real *x) { x[0] = 0; x[1] = 0; }
ZCINLINE void rv2_copy(real *x, const real *src) { x[0] = src[0]; x[1] = src[1]; }

ZCINLINE real rv2_sqr(const real *x) { return x[0]*x[0]+x[1]*x[1]; }
ZCINLINE real rv2_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]); }

ZCINLINE real *rv2_normalize(real *x)
{
  real r = rv2_norm(x);
  if (r > 0.0) {
    r = 1.0/r;
    x[0] *= r;
    x[1] *= r;
  }
  return x;
}

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
ZCINLINE real *rv2_sinc(real *x, real *dx, real s)
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

/* for in-place difference use rv3_dec */
ZCINLINE real *rv2_diff(real *diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv2_sum2(real *sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv2_nsum2(real *sum, const real *a, const real *b)
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
#endif /* RV2_H__ */

