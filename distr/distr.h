#include "def.h"
#include "util.h"
#ifndef DISTR_H__
#define DISTR_H__

typedef struct {
  double s;
  double sf;  /* mean force */
  double sf2; /* mean force^2 */
  double sfr; /* reference mean force */
  double sdv;
} distrsum_t;

typedef struct {
  int xcnt;
  double xmin, xmax, dx;
  distrsum_t *arr;
} distr_t;

INLINE distr_t *distr_open(double xmin, double xmax, double dx)
{
  distr_t *d;
  xnew(d, 1);
  d->xmin = xmin;
  d->xmax = xmax;
  d->dx = dx;
  d->xcnt = (int)((xmax - xmin)/dx + 0.999999f);
  xnew(d->arr, d->xcnt + 1);
  return d;
}

INLINE void distr_close(distr_t *d)
{
  if (!d) return;
  free(d->arr);
  free(d);
}

INLINE void distr_add(distr_t *d, double x, double f, double fr, double dv, double w)
{
  if (x >= d->xmin && x <= d->xmax) {
    distrsum_t *ds = d->arr + (int)((x - d->xmin)/d->dx);
    ds->s += w;
    ds->sf += w * f;
    ds->sf2 += w * f * f;
    ds->sfr += w * fr;
    ds->sdv += w * dv;
  }
}

INLINE int distr_save(distr_t *d, const char *fn)
{
  FILE *fp;
  int i, i0, i1, n = d->xcnt;
  distrsum_t *ds = d->arr;

  for (i0 = 0; i0 < n; i0++)
    if (fabs(ds[i0].s) > 0.) break;
  for (i1 = n - 1; i1 >= i0; i1--)
    if (fabs(ds[i1].s) > 0.) break;

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %.14e %.14e\n", d->xcnt, d->xmin, d->dx);

  for (i = i0; i <= i1; i++) {
    double f = 0.0, f2 = 0.0, fr = 0.0, dv = 0.0;
    ds = d->arr + i;
    if (ds->s > 0.) {
      f = ds->sf / ds->s;
      f2 = ds->sf2 / ds->s;
      fr = ds->sfr / ds->s;
      dv = ds->sdv / ds->s;
    }
    fprintf(fp, "%g %.14e %.14e %.14e %.14e %.14e\n",
      d->xmin + (i + .5) * d->dx, ds->s, f, f2, fr, dv);
  }
  fclose(fp);
  return 0;
}

#endif
