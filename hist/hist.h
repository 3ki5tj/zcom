#include "util.h"
#ifndef HIST_H__
#define HIST_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define HIST_VERBOSE    0x0001

#define HIST_ADDAHALF   0x0010
#define HIST_NOZEROES   0x0020
#define HIST_KEEPLEFT   0x0040
#define HIST_KEEPRIGHT  0x0080
#define HIST_KEEPLEFT2  0x0040
#define HIST_KEEPRIGHT2 0x0080
#define HIST_KEEPEDGE   (HIST_KEEPLEFT|HIST_KEEPRIGHT|HIST_KEEPLEFT2|HIST_KEEPRIGHT2)
#define HIST_KEEPHIST   0x0100
#define HIST_OVERALL    0x0200

#define HIST_ADDITION   0x1000

/* old names */
#define wdist(h,m,n,x0,dx,fn) wdistex(h,m,n,x0,dx,HIST_ADDAHALF|HIST_VERBOSE,fn)
#define wdistex histsave

#define histsave(h,rows,n,xmin,dx,flags,fname) \
  histsavex((const double *)h,rows,n,xmin,dx,flags,NULL,NULL,NULL,fname)

INLINE int histsavex(const double *h, int rows, int n, double xmin, double dx,
    unsigned flags, int (*fwheader)(FILE *, void *),
    double (*fnorm)(int, int, double, double, void *),
    void *pdata, const char *fname);

#define histload(h,rows,n,xmin,dx,flags,fname) \
  histloadx((double *)h,rows,n,xmin,dx,flags,NULL,NULL,NULL,fname)

INLINE int histloadx(double *hist, int rows, int n, double xmin, double dx,
    unsigned flags, int (*frheader)(const char *, void *),
    double (*fnorm)(int, int, double, double, void *),
    void *pdata, const char *fn);

INLINE int histadd(const double *x, double w, double *h, int rows,
    int n, double xmin, double dx, unsigned flags);

/* object oriented wrapper functions */
typedef struct {
  int rows;
  int n;
  double xmin;
  double dx;
  double *arr;
  int (*fwheader)(FILE *, void *);
  int (*frheader)(const char *, void *);
  double (*fnorm)(int, int, double, double, void *);
} hist_t;

typedef hist_t hs_t;

#define hs_open(m, x0, x1, dx) hs_openx(m, x0, x1, dx, NULL, NULL, NULL)
#define hs_open1(x0, x1, dx) hs_open(1, x0, x1, dx)
#define hs_save(hs,fn,flags) hs_savex(hs, fn, NULL, flags)
#define hs_load(hs,fn,flags) hs_loadx(hs, fn, NULL, flags)
#define hs_clear(hs) dblcleararr(hs->arr, hs->rows * hs->n)

INLINE hist_t *hs_openx(int rows, double xmin, double xmax, double dx,
    int (*fwh)(FILE *, void *), int (*frh)(const char*, void *),
    double (*fnorm)(int, int, double, double, void *))
{
  hist_t *hs;

  xnew(hs, 1);
  hs->rows = rows;
  hs->xmin = xmin;
  hs->dx   = dx;
  hs->n = (int)((xmax - xmin)/dx + 0.99999999);
  xnew(hs->arr, hs->n*hs->rows);
  hs->fwheader = fwh;
  hs->frheader = frh;
  hs->fnorm = fnorm;
  return hs;
}

INLINE void hs_close(hist_t *hs)
{
  if (!hs) return;
  if (hs->arr) free(hs->arr);
  memset(hs, 0, sizeof(*hs));
  free(hs);
}

INLINE void hs_check(const hist_t *hs)
{
  die_if (hs == NULL, "hs is %p", (const void *) hs);
  die_if (hs->arr == NULL || hs->rows == 0 || hs->n == 0,
    "hist: arr %p rows %d n %d\n", (const void *)(hs->arr), hs->rows, hs->n);
}

INLINE int hs_savex(const hist_t *hs, const char *fn, void *pdata, unsigned flags)
{
  hs_check(hs);
  return histsavex(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags,
      hs->fwheader, hs->fnorm, pdata, fn);
}

INLINE int hs_loadx(hist_t *hs, const char *fn, void *pdata, unsigned flags)
{
  hs_check(hs);
  return histloadx(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags,
      hs->frheader, hs->fnorm, pdata, fn);
}

INLINE int hs_add(hist_t *hs, const double *x, double w, unsigned flags)
{
  hs_check(hs);
  return histadd(x, w, hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags);
}

#define hs_add1ez(hs, x, flags) hs_add1(hs, 0, x, 1, flags)

INLINE int hs_add1(hist_t *hs, int r, double x, double w, unsigned flags)
{
  hs_check(hs);
  die_if (r >= hs->rows || r < 0, "bad row index %d\n", r);
  return histadd(&x, w, hs->arr + r*hs->n, 1, hs->n, hs->xmin, hs->dx, flags);
}

/* two dimensional version */
INLINE int hist2save(const double *h, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fname);
INLINE int hist2load(double *hist, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fname);
INLINE int hist2add(const double *xarr, const double *yarr, int skip,
    double w, double *h, int rows,
    int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags);

typedef struct {
  int rows;
  int n, m;
  double xmin, ymin;
  double dx, dy;
  double *arr, *dumptr;
} hist2_t;

typedef hist2_t hs2_t;

#define hs2_clear(hs2) dblcleararr(hs2->arr, hs2->rows * hs2->n * hs2->m)

#define hs2_open1(xmin, xmax, dx, ymin, ymax, dy) \
  hs2_open(1, xmin, xmax, dx, ymin, ymax, dy)
#define hs2_opensqr(rows, xmin, xmax, dx) \
  hs2_open(rows, xmin, xmax, dx, xmin, xmax, dx)
#define hs2_opensqr1(xmin, xmax, dx) \
  hs2_opensqr(1, xmin, xmax, dx)

INLINE hist2_t *hs2_open(int rows, double xmin, double xmax, double dx,
    double ymin, double ymax, double dy)
{
  hist2_t *hs2;

  xnew(hs2, 1);
  hs2->rows = rows;
  hs2->xmin = xmin;
  hs2->dx   = dx;
  hs2->n    = (int)((xmax - xmin)/dx + 0.99999999);
  hs2->ymin = ymin;
  hs2->dy   = dy;
  hs2->m    = (int)((ymax - ymin)/dy + 0.99999999);
  xnew(hs2->arr, hs2->n * hs2->m * hs2->rows);
  return hs2;
}

INLINE void hs2_close(hist2_t *hs2)
{
  if (hs2) {
    if (hs2->arr) free(hs2->arr);
    memset(hs2, 0, sizeof(*hs2));
    free(hs2);
  }
}

INLINE void hs2_check(const hist2_t *hs)
{
  die_if (hs == NULL, "hist2 is %p", (const void *) hs);
  die_if (hs->arr == NULL || hs->rows == 0 || hs->n == 0 || hs->m == 0,
    "hist2: arr %p rows %d n %d m %d\n",
    (const void *)(hs->arr), hs->rows, hs->n, hs->m);
}

INLINE int hs2_save(const hist2_t *hs, const char *fn, unsigned flags)
{
  hs2_check(hs);
  return hist2save(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx,
      hs->m, hs->ymin, hs->dy, flags, fn);
}

INLINE int hs2_load(hist2_t *hs, const char *fn, unsigned flags)
{
  hs2_check(hs);
  return hist2load(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx,
      hs->m, hs->ymin, hs->dy, flags, fn);
}

INLINE int hs2_add(hist2_t *hs, const double *x, const double *y, int skip,
    double w, unsigned flags)
{
  hs2_check(hs);
  return hist2add(x, y, skip, w, hs->arr, hs->rows, hs->n, hs->xmin, hs->dx,
      hs->m, hs->ymin, hs->dy, flags);
}

#define hs2_add1ez(hs, x, y, flags) hs2_add1(hs, 0, x, y, 1.0, flags)

INLINE int hs2_add1(hist2_t *hs, int r, double x, double y,
    double w, unsigned flags)
{
  hs2_check(hs);
  return hist2add(&x, &y, 1, w, hs->arr+r*hs->n*hs->m, 1,
     hs->n, hs->xmin, hs->dx, hs->m, hs->ymin, hs->dy, flags);
}

#endif

