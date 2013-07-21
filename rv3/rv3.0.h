#include "util.h"
#include "rng.h"
#ifndef RV3_H__
#define RV3_H__

#ifndef FV3_T
#define FV3_T fv3_t
typedef float fv3_t[3];
#endif

#ifndef DV3_T
#define DV3_T dv3_t
typedef double dv3_t[3];
#endif

#ifndef RV3_T
#define RV3_T rv3_t
typedef real rv3_t[3];
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>



#define rv3_print(r, nm, fmt, nl) rv3_fprint(stdout, r, nm, fmt, nl)

INLINE void rv3_fprint(FILE *fp, real *r, const char *nm, const char *fmt, int nl)
{
  int i;

  if (nm) fprintf(fp, "%s: ", nm);
  for (i = 0; i < 3; i++) fprintf(fp, fmt, r[i]);
  fprintf(fp, "%c", nl ? '\n' : ';');
}



/* due to possible pointer overlap, 'const' are not add to some parameters */

INLINE real *rv3_make(real *x, real a, real b, real c)
{
  x[0] = a;
  x[1] = b;
  x[2] = c;
  return x;
}



#define rv3_makev(rv, v) rv3_make(rv, (real) v[0], (real) v[1], (real) v[2])



#define rv3_zero(x) rv3_make(x, 0, 0, 0)



INLINE real *rv3_copy(real * RESTRICT x, const real *src)
{
  x[0] = src[0];
  x[1] = src[1];
  x[2] = src[2];
  return x;
}



/* use macro to avoid const qualifier of src */
#define rv3_ncopy(x, src, n) memcpy(x, src, 3*n*sizeof(real))



INLINE void rv3_swap(real * RESTRICT x, real * RESTRICT y)
{
  real z[3];
  rv3_copy(z, x);
  rv3_copy(x, y);
  rv3_copy(y, z);
}



INLINE real rv3_sqr(const real *x)
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}



INLINE real rv3_norm(const real *x)
{
  return (real) sqrt( rv3_sqr(x) );
}



/* if x == y, try to use sqr */
INLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}



INLINE real *rv3_cross(real * RESTRICT z, const real *x, const real *y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}



INLINE real *rv3_neg(real *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  x[2] = -x[2];
  return x;
}



INLINE real *rv3_neg2(real * RESTRICT nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  nx[2] = -x[2];
  return nx;
}



INLINE real *rv3_inc(real * RESTRICT x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  x[2] += dx[2];
  return x;
}



INLINE real *rv3_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  x[2] -= dx[2];
  return x;
}



INLINE real *rv3_sinc(real * RESTRICT x, const real *dx, real s)
{
  x[0] += dx[0] * s;
  x[1] += dx[1] * s;
  x[2] += dx[2] * s;
  return x;
}



INLINE real *rv3_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}




/* if y == x, just use smul */
INLINE real *rv3_smul2(real * RESTRICT y, const real *x, real s)
{
  y[0] = x[0] * s;
  y[1] = x[1] * s;
  y[2] = x[2] * s;
  return y;
}



INLINE real *rv3_normalize(real *x)
{
  real r = rv3_norm(x);
  if (r > 0.f) rv3_smul(x, 1.f/r);
  return x;
}



INLINE real *rv3_makenorm(real *v, real x, real y, real z)
{
  return rv3_normalize( rv3_make(v, x, y, z) );
}



/* for in-place difference use rv3_dec */
INLINE real *rv3_diff(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
  return c;
}



/* distance^2 between a and b */
INLINE real rv3_dist2(const real *a, const real *b)
{
  real d[3];
  return rv3_sqr(rv3_diff(d, a, b));
}



/* distance between a and b */
INLINE real rv3_dist(const real *a, const real *b)
{
  return (real) sqrt(rv3_dist2(a, b));
}



/* c = a + b, for in-place addition use rv3_inc */
INLINE real *rv3_add(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
  return c;
}



/* c = - a - b */
INLINE real *rv3_nadd(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = - a[0] - b[0];
  c[1] = - a[1] - b[1];
  c[2] = - a[2] - b[2];
  return c;
}



/* c = a + b * s */
INLINE real *rv3_sadd(real * RESTRICT c, const real *a, const real *b, real s)
{
  c[0] = a[0] + b[0] * s;
  c[1] = a[1] + b[1] * s;
  c[2] = a[2] + b[2] * s;
  return c;
}



/* c = a * s1 + b * s2 */
INLINE real *rv3_lincomb2(real * RESTRICT c, const real *a, const real *b,
    real s1, real s2)
{
  c[0] = a[0] * s1 + b[0] * s2;
  c[1] = a[1] * s1 + b[1] * s2;
  c[2] = a[2] * s1 + b[2] * s2;
  return c;
}



/* angle and gradients of cos(x1-x2-x3) */
INLINE real rv3_cosang(const real *x1, const real *x2, const real *x3,
    real * RESTRICT g1, real * RESTRICT g2, real * RESTRICT g3)
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
INLINE real rv3_ang(const real *x1, const real *x2, const real *x3,
    real * RESTRICT g1, real * RESTRICT g2, real * RESTRICT g3)
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



/* vertical distance from `x' to the line extended by `a' and `b' */
INLINE real rv3_vdist(const real *x, const real *a, const real *b)
{
  real nm[3], d[3], dot;

  rv3_diff(d, x, a);
  rv3_normalize(rv3_diff(nm, a, b));
  dot = rv3_dot(d, nm);
  return rv3_norm(rv3_sinc(d, nm, -dot));
}



/* signed distance from x to the plane extended by a, b, c */
INLINE real rv3_vpdist(const real *x, const real *a, const real *b, const real *c)
{
  real u[3], v[3], m[3];

  rv3_diff(u, b, a);
  rv3_diff(v, c, b);
  rv3_normalize(rv3_cross(m, u, v));
  rv3_diff(u, x, a);
  return rv3_dot(u, m);
}



/* light weight dihedral */
INLINE real rv3_dih(const real *xi, const real *xj, const real *xk, const real *xl,
    real * RESTRICT gi, real * RESTRICT gj, real * RESTRICT gk, real * RESTRICT gl)
{
  real tol, phi, cosphi = 1.f;
  real nxkj, nxkj2, m2, n2;
  real xij[3], xkj[3], xkl[3], uvec[3], vvec[3], svec[3];
  real m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */

  rv3_diff(xij, xi, xj);
  rv3_diff(xkj, xk, xj);
  rv3_diff(xkl, xk, xl);
  nxkj2 = rv3_sqr(xkj);
  nxkj = (real) sqrt(nxkj2);
  tol = (sizeof(real) == sizeof(float)) ? nxkj2 * 6e-8f : nxkj2 * 1e-16f;

  rv3_cross(m, xij, xkj);
  m2 = rv3_sqr(m);
  rv3_cross(n, xkj, xkl);
  n2 = rv3_sqr(n);
  if (m2 > tol && n2 > tol) {
    cosphi = rv3_dot(m, n);
    cosphi /= (real) sqrt(m2 * n2);
    if (cosphi >= 1.f) cosphi = 1.f;
    else if (cosphi < -1.f) cosphi = -1.f;
  }
  phi = (real) acos(cosphi);
  if (rv3_dot(n, xij) < 0.0) phi = -phi;

  /* optionally calculate the gradient */
  if (gi != NULL) {
    if (m2 > tol && n2 > tol) {
      rv3_smul2(gi, m, nxkj/m2);
      rv3_smul2(gl, n, -nxkj/n2);
      rv3_smul2(uvec, gi, rv3_dot(xij, xkj)/nxkj2);
      rv3_smul2(vvec, gl, rv3_dot(xkl, xkj)/nxkj2);
      rv3_diff(svec, uvec, vvec);
      rv3_diff(gj, svec, gi);
      rv3_nadd(gk, svec, gl);
    } else { /* clear the gradients */
      rv3_zero(gi);
      rv3_zero(gj);
      rv3_zero(gk);
      rv3_zero(gl);
    }
  }
  return phi;
}



#define rv3_rnd0(v) rv3_rnd(v, 0, 1)

/* uniformly distributed random vector [a, b) */
INLINE real *rv3_rnd(real *v, real a, real b)
{
  b -= a;
  v[0] = a + b * (real) rnd0();
  v[1] = a + b * (real) rnd0();
  v[2] = a + b * (real) rnd0();
  return v;
}



/* displace `x0' by a random vector in [-a, a)^3 */
INLINE real *rv3_rnddisp(real * RESTRICT x, const real *x0, real a)
{
  x[0] = x0[0] + (real) rnd(-a, a);
  x[1] = x0[1] + (real) rnd(-a, a);
  x[2] = x0[2] + (real) rnd(-a, a);
  return x;
}



/* normally distributed random vector */
INLINE real *rv3_grand0(real *v)
{
  v[0] = (real) grand0();
  v[1] = (real) grand0();
  v[2] = (real) grand0();
  return v;
}



/* normally distributed random vector */
INLINE real *rv3_grand(real *v, real c, real r)
{
  v[0] = c + r * (real) grand0();
  v[1] = c + r * (real) grand0();
  v[2] = c + r * (real) grand0();
  return v;
}



/* displace `x0' by a normally-distributed random vector */
INLINE real *rv3_granddisp(real * RESTRICT x, const real *x0, real a)
{
  x[0] = x0[0] + (real) grand0() * a;
  x[1] = x0[1] + (real) grand0() * a;
  x[2] = x0[2] + (real) grand0() * a;
  return x;
}



/* randomly oriented vector on the sphere of radius r */
#define rv3_rnddir(v, r) rv3_smul(rv3_rnddir0(v), r)

/* randomly oriented vector on the unit sphere */
INLINE real *rv3_rnddir0(real *v)
{
  double a, b, sq, s;

  do { /* projection on the x-y plane */
    a = 2 * rnd0() - 1;
    b = 2 * rnd0() - 1;
    sq = a * a + b * b;
  } while (sq >= 1); /* avoid sin() and cos() */

  s = 2. * sqrt(1 - sq);
  return rv3_make(v, (real) (a * s), (real) (b * s), (real) (1 - 2 * sq));
}




/* randomly orientied vector within the sphere of radius `r' */
#define rv3_rndball(v, r) rv3_smul(rv3_rndball0(v), r)

/* randomly orientied vector within the unit sphere */
INLINE real *rv3_rndball0(real *v)
{
  do {
    rv3_rnd(v, -1, 1);
  } while (rv3_sqr(v) >= 1);
  return v;
}



#endif /* RV3_H__ */

