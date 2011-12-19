#include "def.h"
#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef RC_H__
#define RC_H__

typedef struct { real re, im; } rcomplex_t;

/* complex zero */
rcomplex_t rc_zero = {0.f, 0.f};

/* make complex */
INLINE rcomplex_t rc_make(real re, real im)
{
  rcomplex_t x;
  x.re = re; x.im = im;
  return x;
}

/* complex conjugate */
INLINE rcomplex_t rc_star(rcomplex_t x)
{
  x.im = - x.im;
  return x;
}

/* complex times real */
INLINE rcomplex_t rc_smul(rcomplex_t x, real s)
{
  x.re *= s; x.im *= s;
  return x;
}

/* complex divided by real  */
INLINE rcomplex_t rc_sdiv(rcomplex_t x, real s)
{
  x.re /= s; x.im /= s;
  return x;
}

/* complex addition */
INLINE rcomplex_t rc_add(rcomplex_t x, rcomplex_t y)
{
  x.re += y.re; x.im += y.im;
  return x;
}

/* complex multiplication */
INLINE rcomplex_t rc_mul(rcomplex_t x, rcomplex_t y)
{
  rcomplex_t z;
  z.re = x.re * y.re - x.im * y.im;
  z.im = x.re * y.im + x.im * y.re;
  return z;
}

#endif

