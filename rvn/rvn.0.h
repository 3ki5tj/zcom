#include "util.h"
#include "rng.h"
#include "rv2.h"
#include "rv3.h"
#ifndef RVN_H__
#define RVN_H__
/* D-dimensional real vector */



#ifndef D
#define D 3
#endif


#ifndef DM /* aligned array dimension */
#ifdef RVN_ALIGN
/* double arrays of multiples of 4 can use -xhost options in Intel C
 * and avoid possible compiler optimization mistakes */
#define DM ((D + (RVN_ALIGN) - 1) / (RVN_ALIGN) * (RVN_ALIGN))
#else
/* by default, we set D == DM */
#define DM D
#endif
#endif



#ifndef FVN_T
#define FVN_T fvn_t
typedef float fvn_t[DM];
#endif

#ifndef DVN_T
#define DVN_T dvn_t
typedef double dvn_t[DM];
#endif

#ifndef RVN_T
#define RVN_T rvn_t
typedef real rvn_t[DM];
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
    x[i] = (real) va_arg(vl, double);
#if D < DM
  for (i = D; i < DM; i++) x[i] = (real) 0;
#endif
  va_end(vl);
  return x;
}



INLINE real *rvn_zero(real *x)
{
  int i;

  /* clear all the way up to DM */
  for (i = 0; i < DM; i++) x[i] = 0;
  return x;
}



INLINE real *rvn_copy(real *x, const real *src)
{
  int i;

  for (i = 0; i < DM; i++)
    x[i] = src[i];
  return x;
}



/* use macro to avoid const qualifier of src */
#define rvn_ncopy(x, src, n) memcpy(x, src, D*n*sizeof(real))



INLINE void rvn_swap(real * RESTRICT x, real * RESTRICT y)
{
  rvn_t z;
  rvn_copy(z, x);
  rvn_copy(x, y);
  rvn_copy(y, z);
}


#define rvn_sqr(x) rvn_dot(x, x)
#define rvn_norm(x) (real) sqrt(rvn_sqr(x))

#if DM == 2
#define rvn_dot(x, y) rv2_dot(x, y)
#elif DM == 3
#define rvn_dot(x, y) rv3_dot(x, y)
#else
INLINE real rvn_dot(const real *x, const real *y)
{
  int i;
  real dot = 0;

  /* assuming x[i], y[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) dot += x[i] * y[i];
  return dot;
}
#endif



INLINE real *rvn_neg(real * RESTRICT x)
{
  int i;

  /* assuming x[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) x[i] = -x[i];
  return x;
}



INLINE real *rvn_neg2(real * RESTRICT nx, const real *x)
{
  int i;

  /* assuming x[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) nx[i] = -x[i];
  return nx;
}



INLINE real *rvn_inc(real * RESTRICT x, const real *dx)
{
  int i;

  /* assuming dx[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) x[i] += dx[i];
  return x;
}



INLINE real *rvn_dec(real * RESTRICT x, const real *dx)
{
  int i;

  /* assuming dx[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) x[i] -= dx[i];
  return x;
}



INLINE real *rvn_sinc(real * RESTRICT x, const real *dx, real s)
{
  int i;

  /* assuming dx[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) x[i] += s * dx[i];
  return x;
}



INLINE real *rvn_smul(real *x, real s)
{
  int i;

  for (i = 0; i < DM; i++) x[i] *= s;
  return x;
}



INLINE real *rvn_smul2(real * RESTRICT y, const real *x, real s)
{
  int i;

  /* assuming x[i] == 0 for i >= D */
  for (i = 0; i < DM; i++) y[i] = x[i] * s;
  return y;
}



INLINE real *rvn_normalize(real *x)
{
  real r = rvn_norm(x);
  if (r > (real) 0) rvn_smul(x, (real) 1./r);
  return x;
}



INLINE real *rvn_makenorm(real *v, ...)
{
  int i;
  va_list vl;

  va_start(vl, v);
  for (i = 0; i < D; i++)
    v[i] = (real) va_arg(vl, double);
  va_end(vl);
#if D < DM
  for (i = D; i < DM; i++) v[i] = (real) 0;
#endif
  return rvn_normalize( v );
}



/* for in-place difference use rvn_dec */
INLINE real *rvn_diff(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  /* assuming a[i] = b[i] = 0 for i >= D */
  for (i = 0; i < DM; i++) c[i] = a[i] - b[i];
  return c;
}



/* distance^2 between a and b */
#if DM == 2
#define rvn_dist2(a, b) rv2_dist2(a, b)
#elif DM == 3
#define rvn_dist2(a, b) rv3_dist2(a, b)
#else
INLINE real rvn_dist2(const real *a, const real *b)
{
  rvn_t d;
  return rvn_sqr(rvn_diff(d, a, b));
}
#endif



/* distance between a and b */
INLINE real rvn_dist(const real *a, const real *b)
{
  return (real) sqrt(rvn_dist2(a, b));
}



/* c = a + b, for in-place addition use rvn_inc */
INLINE real *rvn_add(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  /* assuming a[i] = b[i] = 0 for i >= D */
  for (i = 0; i < DM; i++) c[i] = a[i] + b[i];
  return c;
}



/* c = - a - b */
INLINE real *rvn_nadd(real * RESTRICT c, const real *a, const real *b)
{
  int i;

  /* assuming a[i] = b[i] = 0 for i >= D */
  for (i = 0; i < DM; i++) c[i] = - a[i] - b[i];
  return c;
}



/* c = a + b * s */
INLINE real *rvn_sadd(real * RESTRICT c, const real *a, const real *b, real s)
{
  int i;

  /* assuming a[i] = b[i] = 0 for i >= D */
  for (i = 0; i < DM; i++) c[i] = a[i] + b[i] * s;
  return c;
}



/* c = a * s1 + b * s2 */
INLINE real *rvn_lincomb2(real * RESTRICT c, const real *a, const real *b, real s1, real s2)
{
  int i;

  /* assuming a[i] = b[i] = 0 for i >= D */
  for (i = 0; i < DM; i++) c[i] = a[i] * s1 + b[i] * s2;
  return c;
}



/* cosine of the angle of x1-x2-x3 */
INLINE real rvn_cosang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  rvn_t a, b;
  real ra, rb, dot;

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
  rvn_t nm, d;
  real dot;

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
#if D < DM
  for (i = D; i < DM; i++) v[i] = (real) 0;
#endif
  return v;
}



#if DM == 2
#define rvn_rnddisp(x, x0, a) rv2_rnddisp(x, x0, a)
#elif DM == 3
#define rvn_rnddisp(x, x0, a) rv3_rnddisp(x, x0, a)
#else
/* displace `x0' by a random vector in [-a, a)^D */
INLINE real *rvn_rnddisp(real * RESTRICT x, const real *x0, real a)
{
  int i;

  for (i = 0; i < D; i++)
    x[i] = x0[i] + (real) rnd(-a, a);
#if D < DM
  for (i = D; i < DM; i++) x[i] = (real) 0;
#endif
  return x;
}
#endif



/* normally distributed random vector */
INLINE real *rvn_grand0(real *v)
{
  int i;

  for (i = 0; i < D; i++)
    v[i] = (real) grand0();
#if D < DM
  for (i = D; i < DM; i++) v[i] = (real) 0;
#endif
  return v;
}



/* normally distributed random vector */
INLINE real *rvn_grand(real *v, real c, real r)
{
  int i;

  for (i = 0; i < D; i++)
    v[i] = c + r * (real) grand0();
#if D < DM
  for (i = D; i < DM; i++) v[i] = (real) 0;
#endif
  return v;
}



/* displace `x0' by a normally-distributed random vector */
#if DM == 2
#define rvn_granddisp(x, x0, a) rv2_granddisp(x, x0, a)
#elif DM == 3
#define rvn_granddisp(x, x0, a) rv3_granddisp(x, x0, a)
#else
INLINE real *rvn_granddisp(real * RESTRICT x, const real *x0, real a)
{
  int i;

  for (i = 0; i < D; i++)
    x[i] = x0[i] + (real) grand0() * a;
#if D < DM
  for (i = D; i < DM; i++) x[i] = (real) 0;
#endif
  return x;
}
#endif



/* randomly oriented vector on the sphere of radius r */
#define rvn_rnddir(v, r) rvn_smul(rvn_rnddir0(v), r)

/* randomly oriented vector on the unit sphere */
INLINE real *rvn_rnddir0(real *v)
{
#if D < 5
  while ( rvn_sqr(rvn_rnd(v, -1, 1)) >= 1 ) ;
  return rvn_normalize(v);
#else
  /* if D >= 5, normal distribution is faster */
  return rvn_normalize( rvn_grand0(v) );
#endif
}



/* randomly oriented vector within the sphere of radius `r' */
#define rvn_rndball(v, r) rvn_smul(rvn_rndball0(v), r)

/* randomly vector within the unit sphere */
INLINE real *rvn_rndball0(real *v)
{
#if D < 5
  while ( rvn_sqr( rvn_rnd(v, -1, 1) ) >= 1 ) ;
  return v;
#else
  real r = (real) pow(rnd0(), 1.0/D), nm;
  /* first obtain a orientation */
  while ( (nm = rvn_norm(rvn_grand0(v))) <= 1e-8 ) ;
  /* the probability density rho(r) ~ r^(D - 1), so the cumulative
   * distribution P(r) = r^D, and r is obtained from P(r)^(1/D) */
  return rvn_smul(v, r/nm);
#endif
}



#endif /* RVN_H__ */

