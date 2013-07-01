#include "util.h"
#include "rng.h"
#ifndef RV2_H__
#define RV2_H__
#include <stdio.h>
#include <string.h>
#include <math.h>



#ifndef FV2_T
#define FV2_T fv2_t
typedef float fv2_t[2];
#endif

#ifndef DV2_T
#define DV2_T dv2_t
typedef double dv2_t[2];
#endif

#ifndef RV2_T
#define RV2_T rv2_t
typedef real rv2_t[2];
#endif



#define rv2_print(r, nm, fmt, nl) rv2_fprint(stdout, r, nm, fmt, nl)

INLINE void rv2_fprint(FILE *fp, const real *r, const char *nm,
    const char *fmt, int nl)
{
  int i;
  if (nm) fprintf(fp, "%s: ", nm);
  for (i = 0; i < 2; i++)
    fprintf(fp, fmt, r[i], nl);
  fprintf(fp, "%c", (nl ? '\n' : ';'));
}




/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */
INLINE real *rv2_make(real *x, real a, real b)
{
  x[0] = a;
  x[1] = b;
  return x;
}



#define rv2_makev(rx, x) rv2_make(rx, (real) x[0], (real) x[1])



#define rv2_zero(x) rv2_make(x, 0, 0)



INLINE real *rv2_copy(real *x, const real *src)
{
  x[0] = src[0];
  x[1] = src[1];
  return x;
}



/* use macro to avoid const qualifier of src */
#define rv2_ncopy(x, src, n) memcpy(x, src, 2*n*sizeof(real))



INLINE void rv2_swap(real * RESTRICT x, real * RESTRICT y)
{
  real z[2];
  rv2_copy(z, x);
  rv2_copy(x, y);
  rv2_copy(y, z);
}



INLINE real rv2_sqr(const real *x)
{
  return x[0] * x[0] + x[1] * x[1];
}


INLINE real rv2_norm(const real *x)
{
  return (real) sqrt(rv2_sqr(x));
}



INLINE real rv2_dot(const real *x, const real *y)
{
  return x[0] * y[0] + x[1] * y[1];
}



INLINE real rv2_cross(const real *x, const real *y)
{
  return x[0] * y[1] - x[1] * y[0];
}



INLINE real *rv2_neg(real * RESTRICT x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  return x;
}



INLINE real *rv2_neg2(real * RESTRICT nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  return nx;
}



INLINE real *rv2_inc(real * RESTRICT x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  return x;
}



INLINE real *rv2_dec(real * RESTRICT x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  return x;
}



INLINE real *rv2_sinc(real * RESTRICT x, const real *dx, real s)
{
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  return x;
}



INLINE real *rv2_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  return x;
}



INLINE real *rv2_smul2(real * RESTRICT y, const real *x, real s)
{
  y[0] = x[0] * s;
  y[1] = x[1] * s;
  return y;
}



INLINE real *rv2_normalize(real *x)
{
  real r = rv2_norm(x);
  if (r > 0.f) rv2_smul(x, 1.f/r);
  return x;
}



INLINE real *rv2_makenorm(real *v, real x, real y)
{
  return rv2_normalize( rv2_make(v, x, y) );
}



/* for in-place difference use rv2_dec */
INLINE real *rv2_diff(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  return c;
}

/* distance^2 between a and b */
INLINE real rv2_dist2(const real *a, const real *b)
{
  real d[2];
  return rv2_sqr(rv2_diff(d, a, b));
}



/* distance between a and b */
INLINE real rv2_dist(const real *a, const real *b)
{
  return (real) sqrt(rv2_dist2(a, b));
}



/* c = a + b, for in-place addition use rv2_inc */
INLINE real *rv2_add(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  return c;
}



/* c = - a - b */
INLINE real *rv2_nadd(real * RESTRICT c, const real *a, const real *b)
{
  c[0] = - a[0] - b[0];
  c[1] = - a[1] - b[1];
  return c;
}



/* c = a + b * s */
INLINE real *rv2_sadd(real * RESTRICT c, const real *a, const real *b, real s)
{
  c[0] = a[0] + b[0] * s;
  c[1] = a[1] + b[1] * s;
  return c;
}



/* c = a * s1 + b * s2 */
INLINE real *rv2_lincomb2(real * RESTRICT c, const real *a, const real *b, real s1, real s2)
{
  c[0] = a[0] * s1 + b[0] * s2;
  c[1] = a[1] * s1 + b[1] * s2;
  return c;
}



/* cosine of the angle of x1-x2-x3 */
INLINE real rv2_cosang(const real *x1, const real *x2, const real *x3,
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
INLINE real rv2_ang(const real *x1, const real *x2, const real *x3,
    real * RESTRICT g1, real * RESTRICT g2, real * RESTRICT g3)
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
INLINE real rv2_vdist(const real *x, const real *a, const real *b)
{
  real nm[2], d[2], dot;

  rv2_diff(d, x, a);
  rv2_normalize(rv2_diff(nm, a, b));
  dot = rv2_dot(d, nm);
  return rv2_norm(rv2_sinc(d, nm, -dot));
}



/* determinant of a 2x2 matrix */
INLINE real rm2_det(real a[2][2])
{
  return a[0][0]*a[1][1] - a[0][1]*a[1][0];
}



/* inverse matrix b = a^(-1) */
INLINE void rm2_inv(real b[2][2], real a[2][2])
{
  real det = rm2_det(a);
  if (fabs(det) < 1e-30) det = (det < 0) ? -1e-30f: 1e-30f;
  b[0][0] =  a[1][1]/det;
  b[0][1] = -a[0][1]/det;
  b[1][0] = -a[1][0]/det;
  b[1][1] =  a[0][0]/det;
}



#define rv2_rnd0() rv2_rnd(v, 0, 1)

/* uniformly distributed random vector [a, a + b) */
INLINE real *rv2_rnd(real *v, real a, real b)
{
  b -= a;
  v[0] = a + b * (real) rnd0();
  v[1] = a + b * (real) rnd0();
  return v;
}



/* displace `x0' by a random vector in [-a, a)^2 */
INLINE real *rv2_rnddisp(real * RESTRICT x, const real *x0, real a)
{
  x[0] = x0[0] + (real) rnd(-a, a);
  x[1] = x0[1] + (real) rnd(-a, a);
  return x;
}



/* normally distributed random vector */
#define rv2_grand0(v) rv2_grand(v, 0, 1)
INLINE real *rv2_grand(real *v, real c, real r)
{
  v[0] = c + r * (real) grand0();
  v[1] = c + r * (real) grand0();
  return v;
}



/* randomly oriented vector on the sphere of radius r */
#define rv2_rnddir(v, r) rv2_smul(rv2_rnddir0(v), r)

/* randomly oriented vector on the unit sphere */
INLINE real *rv2_rnddir0(real *v)
{
  return rv2_normalize( rv2_rnd(v, -1, 1) );
}




#define rv2_rndball rv2_rnddisk

#define rv2_rndball0 rv2_rnddisk0

/* randomly orientied vector within the sphere of radius `r' */
#define rv2_rnddisk(v, r) rv2_smul(rv2_rnddisk0(v), r)

/* randomly vector within the unit sphere */
INLINE real *rv2_rnddisk0(real *v)
{
  do {
    rv2_rnd(v, -1, 1);
  } while (rv2_sqr(v) >= 1);
  return v;
}



#endif /* RV2_H__ */

