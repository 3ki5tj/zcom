#include "util.h"
#ifndef RC_H__
#define RC_H__
/* routines for complex numbers
 * the functions are similar to the rv2_ ones
 * but here the structures instead of the pointers are passed to the functions */



typedef struct { float re, im; } fcomplex_t;

typedef struct { double re, im; } dcomplex_t;

typedef struct { real re, im; } rcomplex_t;

#if 0
/* complex zero */
static rcomplex_t rc_zero = {0.f, 0.f};
#endif



/* make complex */
INLINE rcomplex_t rc_make(real re, real im)
{
  rcomplex_t a;
  a.re = re;
  a.im = im;
  return a;
}



/* complex conjugate */
INLINE rcomplex_t rc_conj(rcomplex_t a)
{
  a.im = - a.im;
  return a;
}


/* norm */
#define rc_norm(a) (real) sqrt( rc_norm2(a) )

/* square of the norm */
INLINE real rc_norm2(rcomplex_t a)
{
  return a.re * a.re + a.im * a.im;
}



/* complex times real */
INLINE rcomplex_t rc_smul(rcomplex_t a, real s)
{
  a.re *= s;
  a.im *= s;
  return a;
}



/* complex divided by real  */
INLINE rcomplex_t rc_sdiv(rcomplex_t a, real s)
{
  a.re /= s;
  a.im /= s;
  return a;
}



/* addition */
INLINE rcomplex_t rc_add(rcomplex_t a, rcomplex_t b)
{
  a.re += b.re;
  a.im += b.im;
  return a;
}



/* a + s * b */
INLINE rcomplex_t rc_sadd(rcomplex_t a, rcomplex_t b, real s)
{
  a.re += b.re * s;
  a.im += b.re * s;
  return a;
}



/* multiplication */
INLINE rcomplex_t rc_mul(rcomplex_t a, rcomplex_t b)
{
  rcomplex_t c;
  c.re = a.re * b.re - a.im * b.im;
  c.im = a.re * b.im + a.im * b.re;
  return c;
}



/* a / b */
INLINE rcomplex_t rc_div(rcomplex_t a, rcomplex_t b)
{
  real nm2 = rc_norm2(b);
  return rc_sdiv( rc_mul(a, rc_conj(b)), nm2 );
}

#endif

