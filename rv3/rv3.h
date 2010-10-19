#define ZCINLINE __inline static
#define ZCRESTRICT __restrict 

#ifndef RV3_H__
#define RV3_H__

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

ZCINLINE void rv3_zero(real *x) { x[0] = 0.0f; x[1] = 0.0f; x[2] = 0.0f; }
ZCINLINE void rv3_copy(real *x, const real *src) { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; }

ZCINLINE real rv3_sqr (const real *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
ZCINLINE real rv3_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

ZCINLINE real *rv3_normalize(real *x)
{
  real r = rv3_norm(x);
  if (r > 0.0) {
    r = 1.0f/r;
    x[0] *= r;
    x[1] *= r;
    x[2] *= r;
  }
  return x;
}

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
  x[0] -= x[0];
  x[1] -= x[1];
  x[2] -= x[2];
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
ZCINLINE real *rv3_sinc(real *x, real *dx, real s)
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

/* for in-place difference use rv3_dec */
ZCINLINE real *rv3_diff(real * ZCRESTRICT diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  diff[2] = a[2]-b[2];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv3_sum2(real * ZCRESTRICT sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  sum[2] = a[2]+b[2];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv3_nsum2(real *sum, const real *a, const real *b)
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

#endif /* RV3_H__ */

