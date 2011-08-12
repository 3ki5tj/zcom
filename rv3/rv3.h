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

/* for in-place difference use rv3_dec */
ZCINLINE real *rv3_diff(real * ZCRESTRICT diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  diff[2] = a[2]-b[2];
  return diff;
}

/* distance between a and b */
ZCINLINE real rv3_dist(const real *a, const real *b) 
{
  real d[3]; 
  return rv3_norm(rv3_diff(d, a, b));
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

/* determinant of a 3x3 matrix */
ZCINLINE real mat3_det(real a[3][3])
{
  return a[0][0]*a[1][1]*a[2][2]+2.f*a[0][1]*a[0][2]*a[1][2]
    - (a[0][0]*a[1][2]*a[1][2]+a[1][1]*a[0][2]*a[0][2]+a[2][2]*a[0][1]*a[0][1]);
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
    2.f*(a[0][1]*a[0][1]+a[0][2]*a[0][2]+a[1][2]*a[1][2]))/6.f;
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
  double det, tol = 1e-15f;

  vec[2] = 1.f;
  if (fabs(det = a*b - d*d) > tol) { /* use row 0 and 1 */
    vec[0] = (real)((d*f-b*e)/det);
    vec[1] = (real)((e*d-a*f)/det);
    return rv3_normalize(vec);
  }
  c = m[2][2] - val;
  if (fabs(det = a*f - e*d) > tol) { /* row 1 and 2 */
    vec[0] = (real)((d*c-e*f)/det);
    vec[1] = (real)((e*e-a*c)/det);
    return rv3_normalize(vec);
  }
  if ((det = sqrt(a*a+d*d)) > tol) { /* three-row-degenerate */
    vec[0] = (real)(d/det);
    vec[1] = (real)(-a/det);
    vec[2] = 0.f;
  } else {
    vec[0] = 1.f;
    vec[1] = 0.f;
    vec[2] = 0.f;
  }
  return vec;
}

#endif /* RV3_H__ */

