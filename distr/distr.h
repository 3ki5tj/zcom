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

INLINE void distrsum_clear(distrsum_t *x) { x->s = x->sf = x->sf2 = x->sfr = x->sdv = 0.0; }

typedef struct {
  int n;
  double xmin, xmax, dx;
  distrsum_t *arr;
  double *rho, *lnrho, *mf, *his;
} distr_t;

/* integral identity from [i - m, i + m), set d->mf, d->lnrho, ln->rho */
INLINE void distr_ii0(distr_t *d, int m)
{
  int i, j, j0, j1, n = d->n;
  double x, den, tot, dx = d->dx;
  distrsum_t *ds;

  /* construct ln(rho) from the mean force */
  for (i = 0; i < n; i++) {
    ds = d->arr + i;
    d->mf[i] = (ds->s > 0) ? (ds->sf - ds->sfr) / ds->s : 0.0;
    d->lnrho[i + 1] = d->lnrho[i] + d->mf[i] * dx;
  }

  /* compute the total visits */
  for (tot = 0., i = 0; i < n; i++) tot += d->arr[i].s;

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i */
    j0 = intmax(i - m, 0);
    j1 = intmin(i + m, n);
    for (d->his[i] = 0., j = j0; j < j1; j++)
      d->his[i] += d->arr[j].s;

    /* integrate over exp(lnrho[j] - lnrho[i]) */
    for (den = 0., j = j0; j <= j1; j++) {
      x = exp(d->lnrho[j] - d->lnrho[i]);
      den += (j == j0 || j == j1) ? x * .5 : x;
    }
    den *= tot * dx;
    d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
    d->his[i] /= tot * dx * (j1 - j0);
  }
}

/* Adib-Jazynski identity from [i - m, i + m) 
 * linear function phi(x) is assumed */
INLINE void distr_aj(distr_t *d, int m)
{
  int i, j, j0, j1, n = d->n;
  double x, s0, s1, phi, den, tot, dx = d->dx;
  distrsum_t *ds;

  /* compute the total visits */
  for (tot = 0., i = 0; i < n; i++) tot += d->arr[i].s;

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i */
    j0 = intmax(i - m, 0);
    j1 = intmin(i + m, n);
    for (s0 = s1 = 0., j = j0; j < j1; j++) {
      phi = (j >= i) ? (j + .5 - j1) : (j + .5 - j0);
      ds = d->arr + j;
      s0 += ds->s; /* phi' = 1 */
      s1 += (ds->sf - ds->sfr) * phi;
    }
    den = 2.0 * m * tot * dx;
    d->his[i] = s0/den;
    d->rho[i] = (s0 + s1)/den;
  }
}


INLINE distr_t *distr_open(double xmin, double xmax, double dx)
{
  distr_t *d;
  xnew(d, 1);
  d->xmin = xmin;
  d->xmax = xmax;
  d->dx = dx;
  d->n = (int)((xmax - xmin)/dx + 0.999999f);
  xnew(d->arr, d->n + 1);
  xnew(d->rho, d->n + 1);
  xnew(d->lnrho, d->n + 1);
  xnew(d->mf, d->n + 1);
  xnew(d->his, d->n + 1);
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
  int i, i0 = 0, i1 = d->n, n = d->n;
  distrsum_t *ds = d->arr;

  /* determine the range */
  for (i0 = 0; i0 <= n; i0++)
    if (d->his[i0] > 0. || d->arr[i0].s > 0.)
      break;
  for (; i1 >= i0; i1--)
    if (d->his[i1] > 0. || d->arr[i1].s > 0.)
      break;
  if (i0 > i1) {
    printf("histogram is empty\n");
    return 1;
  }

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d %.14e %.14e\n", 0, d->n, d->xmin, d->dx);

  for (i = i0; i <= i1; i++) {
    double f = 0.0, f2 = 0.0, fr = 0.0, dv = 0.0;
    ds = d->arr + i;
    if (i < n && ds->s > 0.) {
      f = ds->sf / ds->s;
      f2 = ds->sf2 / ds->s - f * f; /* convert to variance */
      fr = ds->sfr / ds->s;
      dv = ds->sdv / ds->s;
    }
    fprintf(fp, "%g %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
      d->xmin + (i + .5) * d->dx, ds->s, f, f2, fr, dv, 
      d->rho[i], d->lnrho[i], d->mf[i], d->his[i]);
  }
  fclose(fp);
  return 0;
}

/* load a previous histogram */
int distr_load(distr_t *d, const char *fn)
{
  FILE *fp;
  char s[1024];
  int ver, i, n = d->n, nlin;
  double x, sm, f, f2, fr, dv, xmin = d->xmin, dx = d->dx;
  distrsum_t *ds;
    
  xfopen(fp, fn, "r", return -1);
  
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  if (4 != sscanf(s, " # %d%d%lf%lf", &ver, &i, &x, &sm)
      || i != n || fabs(x - xmin) > 1e-3 || fabs(sm - dx) > 1e-5) {
    fprintf(stderr, "error: bins = %d, %d, xmin = %g, %g; dx = %g, %g\n", 
        i, n, x, xmin, sm, dx);
    fclose(fp);
    return -1;
  }

  /* only read in raw data */
  for (nlin = 2; fgets(s, sizeof s, fp); nlin++) {
    if (6 != sscanf(s, "%lf%lf%lf%lf%lf%lf", &x, &sm, &f, &f2, &fr, &dv)) goto ERR;
    /* locate the bin */
    i = (int)( (x - xmin)/dx );
    if (i < 0 || i >= n) {
      fprintf(stderr, "bad x %g, xmin %g, i %d, n %d\n", x, xmin, i, n);
      goto ERR;
    }
    ds = d->arr + i;
    ds->s = sm;
    ds->sf = f * sm;
    ds->sf2 = (f2 + f*f)*sm;
    ds->sfr = fr * sm;
    ds->sdv = dv * sm;
  }
  return 0;
ERR:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  /* clear everything on error */
  for (i = 0; i <= n; i++) distrsum_clear(d->arr + i);
  return -1;
}

#endif
