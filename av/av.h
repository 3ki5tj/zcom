#ifndef INLINE
#define INLINE __inline static
#endif
#include "util.h"
#ifndef AV_H__
#define AV_H__
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

typedef struct {
  double s, sx, sx2; /* sum, sum x, variance */
} av_t;

INLINE void av_clear(av_t *av) { av->s = av->sx = av->sx2 = 0; }
INLINE double av_getave(const av_t *av) { return (av && av->s > 0) ? av->sx/av->s : 0; }
INLINE double av_getvar(const av_t *av) { return (av && av->s > 0) ? av->sx2/av->s : 0; }
INLINE double av_getdev(const av_t *av) { return (av && av->s > 0) ? sqrt(av_getvar(av)) : 0; }

/* add a new value to av_t with a weight `w' */
INLINE void av_addw(av_t *av, double x, double w)
{
  double s, sx;

  av->s = (s = av->s) + w;
  av->sx = (sx = av->sx) + x*w;
  if (s <= 0.0) return;
  av->sx2 += (x - sx/s)*(x - av->sx/av->s)*w;
}
#define av_add(av, x) av_addw(av, x, 1)

/* update: sX = sX*gam + X */
INLINE void av_gaddw(av_t *av, double x, double w, double ngam)
{
  double s, sx, del, gam = 1.0 - ngam;

  av->s = (s = av->s)*gam + w;
  av->sx = (sx = av->sx)*gam + w*x;
  if (s <= 0.0) return;
  del = x*s - sx;
  av->sx2 = (av->sx2 + w*del*del/(s*av->s))*gam;
}

#define av_gadd(av, x, ngam) av_gaddw(av, x, 1, ngam)


/* average of n quantities */
typedef struct {
  int n;
  double s;
  double *sx; /* sum */
  double *sx2; /* variance */
} avn_t;

/* open average for n quantities */
INLINE avn_t *avn_open(int n)
{
  avn_t *a;
  int i;

  xnew(a, 1);
  die_if (n <= 0, "avn needs at least n %d >= 1\n", n);
  a->n = n;
  xnew(a->sx, n);
  xnew(a->sx2, n);
  for (a->s = 0, i = 0; i < n; i++)
    a->sx[i] = a->sx2[i] = 0;
  return a;
}

INLINE void avn_close(avn_t *a)
{
  free(a->sx);
  free(a->sx2);
  free(a);
}

/* add values to avn_t with a weight `w'
 * must add all values simultaneously, otherwise a->s is messed up */
INLINE void avn_addwv(avn_t *a, const double *x, double w)
{
  int k, n = a->n;
  double s, sx;

  a->s = (s = a->s) + w;
  for (k = 0; k < n; k++) {
    a->sx[k] = (sx = a->sx[k]) + x[k] * w;
    if (s > 0)
      a->sx2[k] += (x[k] - sx/s) * (x[k] - a->sx[k]/a->s)*w;
  }
}

#define avn_addv(a, x) avn_addwv(a, x, 1)

/* add values to avn_t with a weight `w'
 * use argument list */
INLINE void avn_addw(avn_t *a, double w, ...)
{
  int k, n = a->n;
  double s, sx, x;
  va_list va;

  a->s = (s = a->s) + w;
  va_start(va, w);
  for (k = 0; k < n; k++) {
    x = va_arg(va, double);
    a->sx[k] = (sx = a->sx[k]) + x * w;
    if (s > 0)
      a->sx2[k] += (x - sx/s) * (x - a->sx[k]/a->s)*w;
  }
  va_end(va);
}


/* weighted update: sX = sX*gam + X */
INLINE void avn_gaddwv(avn_t *a, const double *x, double w, double ngam)
{
  int k, n = a->n;
  double s, sx, del, gam = 1.0 - ngam;

  a->s = (s = a->s)*gam + w;
  for (k = 0; k < n; k++) {
    a->sx[k] = (sx = a->sx[k]) * gam + w * x[k];
    if (s > 0.0) {
      del = x[k]*s - sx;
      a->sx2[k] = (a->sx2[k] + w*del*del/(s*a->s))*gam;
    }
  }
}

#define avn_gaddv(a, x, ngam) avn_gaddwv(a, x, 1, ngam)

/* weighted update
 * use argument list */
INLINE void avn_gaddw(avn_t *a, double w, double ngam, ...)
{
  int k, n = a->n;
  double s, sx, x, del, gam = 1 - ngam;
  va_list va;

  a->s = (s = a->s)*gam + w;
  va_start(va, ngam);
  for (k = 0; k < n; k++) {
    x = va_arg(va, double);
    a->sx[k] = (sx = a->sx[k]) * gam + w * x ;
    if (s > 0.0) {
      del = x*s - sx;
      a->sx2[k] = (a->sx2[k] + w*del*del/(s*a->s))*gam;
    }
  }
  va_end(va);
}

/* these macros are only available if we have variable arguments macros */
#ifdef HAVEVAM
#define anv_add(a, ...) avn_addw(a, 1.0, ## __VARARGS__)
#define anv_gadd(a, ngam, ...) avn_gaddw(a, 1.0, ngam, ## __VARARGS__)
#endif

INLINE void avn_clear(avn_t *a)
{
  int i;

  for (a->s = 0, i = 0; i < a->n; i++)
    a->sx[i] = a->sx2[i] = 0;
}

/* get average of quantity k */
INLINE double avn_getave(const avn_t *a, int k)
{
  die_if (k < 0 || k >= a->n, "avn index %d out of range %d\n", k, a->n);
  return (a->s > 0) ? a->sx[k]/a->s : 0; 
}

/* get averages of all quantities */
INLINE double *avn_getaven(const avn_t *a, double *val)
{
  int k;

  if (a->s <= 0.) {
    for (k = 0; k < a->n; k++) val[k] = 0;
  } else {
    for (k = 0; k < a->n; k++)
      val[k] = a->sx[k] / a->s;
  }
  return val;
}

/* get variance of quantity k */
INLINE double avn_getvar(const avn_t *a, int k)
{
  die_if (k < 0 || k >= a->n, "avn index %d out of range %d\n", k, a->n);
  return (a->s > 0) ? a->sx2[k]/a->s : 0; 
}

/* get variances of all quantities */
INLINE double *avn_getvarn(const avn_t *a, double *val)
{
  int k;

  if (a->s <= 0.) {
    for (k = 0; k < a->n; k++) val[k] = 0;
  } else {
    for (k = 0; k < a->n; k++)
      val[k] = a->sx2[k] / a->s;
  }
  return val;
}

/* get standard deviation of quantity k */
INLINE double avn_getdev(const avn_t *a, int k) 
{
  die_if (k < 0 || k >= a->n, "avn index %d out of range %d\n", k, a->n);
  return (a->s > 0) ? sqrt(a->sx2[k]/a->s) : 0; 
}

/* get standard deviations of all quantities */
INLINE double *avn_getdevn(const avn_t *a, double *val)
{
  int k;

  if (a->s <= 0.) {
    for (k = 0; k < a->n; k++) val[k] = 0;
  } else {
    for (k = 0; k < a->n; k++)
      val[k] = sqrt(a->sx2[k] / a->s);
  }
  return val;
}


#endif
