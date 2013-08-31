#include "util.h"
#ifndef AV_H__
#define AV_H__
#include <stdio.h>
#include <stdarg.h>
#include <math.h>



typedef struct {
  double s, sx;
} av0_t; /* simplest averager without variance */



INLINE void av0_clear(av0_t *av) { av->s = 1e-20; av->sx = 0; }
INLINE double av0_getave(const av0_t *av) { return av->sx / av->s; }
INLINE void av0_add(av0_t *av, double x) { av->s += 1; av->sx += x; }
INLINE void av0_addw(av0_t *av, double x, double w) { av->s += w; av->sx += x * w; }




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

/* adaptive averaging: sX = sX * gam + X */
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
  double *x; /* buffer, current x value */
  double *sx, *sxb; /* sum */
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
  xnew(a->x, n);
  xnew(a->sx, n);
  xnew(a->sxb, n);
  xnew(a->sx2, n*n);
  a->s = 0;
  for (i = 0; i < n; i++) a->x[i] = a->sx[i] = a->sxb[i] = 0;
  for (i = 0; i < n * n; i++) a->sx2[i] = 0;
  return a;
}

INLINE void avn_close(avn_t *a)
{
  free(a->x);
  free(a->sx);
  free(a->sxb);
  free(a->sx2);
  free(a);
}

/* add values to avn_t with a weight `w'
 * must add all values simultaneously, otherwise a->s is messed up */
INLINE void avn_addwv(avn_t *a, const double *x, double w)
{
  int k, l, n = a->n;
  double s;

  a->s = (s = a->s) + w;
  for (k = 0; k < n; k++) {
    a->sx[k] = (a->sxb[k] = a->sx[k]) + x[k] * w;
  }
  if (s > 0) { /* update variance */
    for (k = 0; k < n; k++)
      for (l = k; l < n; l++)
        a->sx2[k*n + l] += (x[k] - a->sxb[k]/s) * (x[l] - a->sx[l]/a->s) * w;
  }
}

#define avn_addv(a, x) avn_addwv(a, x, 1)

/* add values to avn_t with a weight `w'
 * use argument list */
INLINE void avn_addw(avn_t *a, double w, ...)
{
  int k;
  va_list va;

  va_start(va, w);
  for (k = 0; k < a->n; k++) {
    a->x[k] = va_arg(va, double);
  }
  va_end(va);
  avn_addwv(a, a->x, w);
}


/* weighted update: sX = sX*gam + X */
INLINE void avn_gaddwv(avn_t *a, const double *x, double w, double ngam)
{
  int k, l, n = a->n;
  double s, gam = 1.0 - ngam;

  a->s = (s = a->s)*gam + w;
  for (k = 0; k < n; k++) {
    a->sx[k] = (a->sxb[k] = a->sx[k]) * gam + w * x[k];
  }
  if (s > 0) { /* update variance */
    for (k = 0; k < n; k++)
      for (l = k; l < n; l++)
        a->sx2[k*n + l] = gam * (a->sx2[k*n + l] +
            w*(x[k]*s - a->sxb[k]) * (x[l]*s - a->sxb[l])/(s*a->s));
  }
}

#define avn_gaddv(a, x, ngam) avn_gaddwv(a, x, 1, ngam)

/* weighted update
 * use argument list */
INLINE void avn_gaddw(avn_t *a, double w, double ngam, ...)
{
  int k;
  va_list va;

  va_start(va, ngam);
  for (k = 0; k < a->n; k++)
    a->x[k] = va_arg(va, double);
  va_end(va);
  avn_gaddwv(a, a->x, w, ngam);
}

/* these macros are only available if we have variable arguments macros */
#ifdef HAVEVAM
#define anv_add(a, ...) avn_addw(a, 1.0, ## __VARARGS__)
#define anv_gadd(a, ngam, ...) avn_gaddw(a, 1.0, ngam, ## __VARARGS__)
#endif

INLINE void avn_clear(avn_t *a)
{
  int i;

  a->s = 0;
  for (i = 0; i < a->n; i++)
    a->x[i] = a->sx[i] = a->sxb[i] = 0;
  for (i = 0; i < a->n * a->n; i++)
    a->sx2[i] = 0;
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

/* get cross variance of quantities k and l */
INLINE double avn_getvar(const avn_t *a, int k, int l)
{
  die_if (k < 0 || k >= a->n || l < 0 || l >= a->n,
      "avn index %d, %d out of range %d\n", k, l, a->n);
  if (k > l) intswap(k, l);
  return (a->s > 0) ? a->sx2[k * a->n + l]/a->s : 0;
}

/* get variances of all quantities */
INLINE double *avn_getvarn(const avn_t *a, double *val)
{
  int k, l, n = a->n;

  if (a->s <= 0.) {
    for (k = 0; k < n * n; k++) val[k] = 0;
  } else {
    for (k = 0; k < n; k++) {
      for (l = k; l < n; l++) {
        val[k*n + l] = a->sx2[k*n + l] / a->s;
        if (l > k) val[l*n + k] = val[k*n + l];
      }
    }
  }
  return val;
}

/* get standard deviation of quantity k */
INLINE double avn_getdev(const avn_t *a, int k)
{
  die_if (k < 0 || k >= a->n, "avn index %d out of range %d\n", k, a->n);
  return (a->s > 0) ? sqrt(a->sx2[k * a->n + k]/a->s) : 0;
}

/* get standard deviations of all quantities */
INLINE double *avn_getdevn(const avn_t *a, double *val)
{
  int k, n = a->n;

  if (a->s <= 0.) {
    for (k = 0; k < n; k++) val[k] = 0;
  } else {
    for (k = 0; k < n; k++)
      val[k] = sqrt(a->sx2[k*n + k] / a->s);
  }
  return val;
}

/* get correlation coefficient between quantities k and l */
INLINE double avn_getcor(const avn_t *a, int k, int l)
{
  int n = a->n;
  die_if (k < 0 || k >= n || l < 0 || l >= n,
      "avn index %d, %d out of range %d\n", k, l, n);
  if (k > l) intswap(k, l);
  return (a->s > 0) ? a->sx2[k*n + l] / sqrt(a->sx2[k*n + k] * a->sx2[l*n + l]) : 0;
}

/* get correlation coefficients among all quantities */
INLINE double *avn_getcorn(const avn_t *a, double *val)
{
  int k, l, n = a->n;

  if (a->s <= 0.) {
    for (k = 0; k < n * n; k++) val[k] = 0;
  } else {
    for (k = 0; k < n; k++) {
      for (l = k; l < n; l++) {
        val[k*n + l] = a->sx2[k*n + l] / sqrt(a->sx2[k*n + k] * a->sx2[l*n + l]);
        if (l > k) val[l*n + k] = val[k*n + l];
      }
    }
  }
  return val;
}

#endif
