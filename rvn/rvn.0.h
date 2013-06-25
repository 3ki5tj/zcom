#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef RESTRICT
#define RESTRICT __restrict
#endif
#include "def.h"
#include "rng.h"
#ifndef RVN_H__
#define RVN_H__
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#ifndef D
#define D 3
#endif


#ifndef FVN_T
#define FVN_T fvn_t
typedef float fvn_t[D];
#endif

#ifndef DVN_T
#define DVN_T dvn_t
typedef double dvn_t[D];
#endif

#ifndef RVN_T
#define RVN_T rvn_t
typedef real rvn_t[D];
#endif



#define rvn_print(r, nm, fmt, nl) rvn_fprint(stdout, r, nm, fmt, nl)

INLINE void rvn_fprint(FILE *fp, const real *r, const char *nm,
    const char *fmt, int nl)
{
  int i;
  if (nm) fprintf(fp, "%s: ", nm);
  for (i = 0; i < D; i++)
    fprintf(fp, fmt, r[i], nl);
  fprintf(fp, "%c", (nl ? '\n' : ';'));
}



INLINE real *rvn_make(real *x, ...)
{
  int i;
  va_list vl;

  va_start(vl, x);
  for (i = 0; i < D; i++)
    x[i] = va_arg(vl, real);
  va_end(vl);
  return x;
}



INLINE real *rvn_zero(real *x)
{
  int i;
  
  for (i = 0; i < D; i++) x[i] = 0;
  return x;
}



INLINE real *rvn_copy(real *x, const real *src)
{
  int i;
  
  for (i = 0; i < D; i++) 
    x[i] = src[i];
  return x;
}



/* use macro to avoid const qualifier of src */
#define rvn_ncopy(x, src, n) memcpy(x, src, 2*n*sizeof(real))



INLINE void rvn_swap(real * RESTRICT x, real * RESTRICT y)
{
  real z[D];
  rvn_copy(z, x);
  rvn_copy(x, y);
  rvn_copy(y, z);
}



INLINE real rvn_sqr(const real *x)
{
  int i;
  real dot = 0;

  for (i = 0; i < D; i++) dot += x[i] * x[i];
  return dot;
}


INLINE real rvn_norm(const real *x)
{
  return (real) sqrt(rvn_sqr(x));
}



INLINE real rvn_dot(const real *x, const real *y)
{
  int i;
  real dot = 0;

  for (i = 0; i < D; i++) dot += x[i] * y[i];
  return dot;
}



INLINE real *rvn_neg(real * RESTRICT x)
{
  int i;

  for (i = 0; i < D; i++) x[i] -= x[i];
  return x;
}



INLINE real *rvn_neg2(real * RESTRICT nx, const real *x)
{
  int i;

  for (i = 0; i < D; i++) nx[i] -= x[i];
  return nx;
}



INLINE real *rvn_inc(real * RESTRICT x, const real *dx)
{
  int i;

  for (i = 0; i < D; i++) x[i] += dx[i];
  return x;
}



INLINE real *rvn_dec(real * RESTRICT x, const real *dx)
{
  int i;

  for (i = 0; i < D; i++) x[i] -= dx[i];
  return x;
}



INLINE real *rvn_sinc(real * RESTRICT x, const real *dx, real s)
{
  int i;

  for (i = 0; i < D; i++) x[i] += s * dx[i];
  return x;
}



INLINE real *rvn_smul(real *x, real s)
{
  int i;

  for (i = 0; i < D; i++) x[i] *= s;
  return x;
}



INLINE real *rvn_smul2(real * RESTRICT y, const real *x, real s)
{
  int i;

  for (i = 0; i < D; i++) y[i] = x[i] * s;
  return y;
}



INLINE real *rvn_normalize(real *x)
{
  real r = rvn_norm(x);
  if (r > 0.f) rvn_smul(x, 1.f/r);
  return x;
}



INLINE real *rvn_makenorm(real *v, real x, real y)
{
  return rvn_normalize( rvn_make(v, x, y) );
}



/* for in-place difference use rvn_dec */
INLINE real *rvn_diff(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  for (i = 0; i < D; i++) c[i] = a[i] - b[i];
  return c;
}

/* distance^2 between a and b */
INLINE real rvn_dist2(const real *a, const real *b)
{
  real d[D];
  return rvn_sqr(rvn_diff(d, a, b));
}



/* distance between a and b */
INLINE real rvn_dist(const real *a, const real *b)
{
  return (real) sqrt(rvn_dist2(a, b));
}



/* c = a + b, for in-place addition use rvn_inc */
INLINE real *rvn_add(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  for (i = 0; i < D; i++) c[i] = a[i] + b[i];
  return c;
}



/* c = - a - b */
INLINE real *rvn_nadd(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  for (i = 0; i < D; i++) c[i] = - a[i] - b[i];
  return c;
}



/* c = a + b * s */
INLINE real *rvn_sadd(real * RESTRICT c, const real *a, const real *b, real s)
{
  int i;

  for (i = 0; i < D; i++)
    c[i] = a[i] + b[i] * s;
  return c;
}



/* c = a * s1 + b * s2 */
INLINE real *rvn_lincomb2(real * RESTRICT c, const real *a, const real *b, real s1, real s2)
{
  int i;
  
  for (i = 0; i < D; i++)
    c[i] = a[i] * s1 + b[i] * s2;
  return c;
}



/* cosine of the angle of x1-x2-x3 */
INLINE real rvn_cosang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real a[D], b[D], ra, rb, dot;

  ra = rvn_norm(rvn_diff(a, x1, x2));
  rvn_smul(a, 1.f/ra);
  rb = rvn_norm(rvn_diff(b, x3, x2));
  rvn_smul(b, 1.f/rb);
  dot = rvn_dot(a, b);
  if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
  if (g1) {
    rvn_lincomb2(g1, b, a, 1.f/ra, -dot/ra);
    rvn_lincomb2(g3, a, b, 1.f/rb, -dot/rb);
    rvn_nadd(g2, g1, g3);
  }
  return dot;
}



/* angle and gradients of x1-x2-x3 */
INLINE real rvn_ang(const real *x1, const real *x2, const real *x3,
    real * RESTRICT g1, real * RESTRICT g2, real * RESTRICT g3)
{
  real dot, sn;
  dot = rvn_cosang(x1, x2, x3, g1, g2, g3);
  sn = (real) sqrt(1 - dot*dot);
  if (sn > 1e-7) sn = -1/sn; else sn = 0.;
  if (g1) {
    rvn_smul(g1, sn);
    rvn_smul(g2, sn);
    rvn_smul(g3, sn);
  }
  return (real) acos(dot);
}



/* vertical distance from x to line a-b */
INLINE real rvn_vdist(const real *x, const real *a, const real *b)
{
  real nm[D], d[D], dot;

  rvn_diff(d, x, a);
  rvn_normalize(rvn_diff(nm, a, b));
  dot = rvn_dot(d, nm);
  return rvn_norm(rvn_sinc(d, nm, -dot));
}



#define rvn_rnd0() rvn_rnd(v, 0, 1)

/* uniformly distributed random vector [a, a + b) */
INLINE real *rvn_rnd(real *v, real a, real b)
{
  int i;

  b -= a;
  for (i = 0; i < D; i++)
    v[i] = a + b * (real) rnd0();
  return v;
}


/* normally distributed random vector */
#define rvn_grand0(v) rvn_grand(v, 0, 1)
INLINE real *rvn_grand(real *v, real c, real r)
{
  int i;

  for (i = 0; i < D; i++)
    v[i] = c + r * (real) grand0();
  return v;
}



/* randomly oriented vector on the sphere of radius r */
#define rvn_rnddir(v, r) rvn_smul(rvn_rnddir0(v), r)

/* randomly oriented vector on the unit sphere */
INLINE real *rvn_rnddir0(real *v)
{
  return rvn_normalize( rvn_rnd(v, -1, 1) );
}




/* randomly orientied vector within the sphere of radius `r' */
#define rvn_rndball(v, r) rvn_smul(rvn_rndball0(v), r)

/* randomly vector within the unit sphere */
INLINE real *rvn_rndball0(real *v)
{
  do {
    rvn_rnd(v, -1, 1);
  } while (rvn_sqr(v) >= 1);
  return v;
}



#endif /* RV2_H__ */

