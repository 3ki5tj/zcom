#include "def.h"
#include "util.h"
#ifndef DISTR_H__
#define DISTR_H__

typedef struct {
  double s;
  double sf;  /* mean force */
  double sf2; /* mean force^2 */
  double sdv; /* average value for the 2nd derivative */
} distrsum_t;

INLINE void distrsum_clear(distrsum_t *x) { x->s = x->sf = x->sf2 = x->sdv = 0.0; }

typedef struct {
  int n;
  double xmin, xmax, dx;
  distrsum_t *arr;
  double *rho, *lnrho, *his;
  double *mf, *mdv;
  int *jl, *jr; /* window size */
  double *err;
  double *hsum; /* histogram sum for comparison */
} distr_t;


/* compute mean force from a single bin or a narrowest symmetric window
 * that contains at least one data point */
INLINE void distr_mf0(distr_t *d)
{
  int i, m = 0, n = d->n;
  double s, sf, sdv;

  for (i = 0; i < n; i++) {
    distrsum_t *ds = d->arr + i;
    s = ds->s;
    sf = ds->sf;
    sdv = ds->sdv;
    /* extend to both sides, stop as soon as we have one data point */
    for (m = 1; m < n && s <= 0.; m++) {
      if (i - m >= 0) { /* extend to the left */
        ds = d->arr + i - m;
        s += ds->s;
        sf += ds->sf;
        sdv += ds->sdv;
      }
      if (i + m < n) { /* extend to the right */
        ds = d->arr + i + m;
        s += ds->s;
        sf += ds->sf;
        sdv += ds->sdv;
      }
    }
    d->mf[i] = (s > 0.) ? (sf / s) : 0.;
    d->mdv[i] = (s > 0.) ? (sdv / s) : 0.;
  }
}

#define distr_mfav(d, m) distr_mfwin(d, m, 0)
#define distr_mfii(d, m) distr_mfwin(d, m, 1)

/* compute mean force from a symmetric window
 * if ii, we use the integral identity, otherwise, a plain average
 * half window size starts from m until the window contains at least one data point
 * linear phi(x) is used for the integral identity
 * improves mean force, slightly smoothness of distribution itself */
INLINE void distr_mfwin(distr_t *d, int m, int ii)
{
  int i, j, j0, j1, mm, n = d->n;
  double s0, s1, s, phi, dx = d->dx;

  distr_mf0(d);
  for (i = 0; i < n; i++) {
    /* start with a window size of m, but extend if nothing is in the window */
    for (mm = m; mm < n; mm++) {
      j0 = intmax(0, i - mm);
      j1 = intmin(n, i + mm + 1);
      for (s0 = s1 = s = 0., j = j0; j < j1; j++) {
        phi = (j == i) ? 0 : (j + .5 - (j < i ? j0 : j1)) * dx;
        s0 += d->mf[j];
        s1 += d->mdv[j] * phi;
        s += 1;
      }
      if (s > 0) break;
    }
    if (ii) s0 += s1; /* add correction second-order derivative */
    d->mf[i] = (s > 0) ? s0/s : 0.0;
  }
}

/* Adib-Jazynski identity from [i - m, i + m)
 * linear function phi(x) is assumed */
INLINE void distr_aj(distr_t *d, int m)
{
  int i, j, j0, j1, n = d->n;
  double s0, s1, phi, den, tot, dx = d->dx;
  distrsum_t *ds;

  /* compute the total visits */
  for (tot = 0., i = 0; i < n; i++) tot += d->arr[i].s;

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i */
    j0 = intmax(i - m, 0);
    j1 = intmin(i + m, n);
    for (s0 = s1 = 0., j = j0; j < j1; j++) {
      phi = (j + .5 - (j < i ? j0 : j1)) * dx;
      ds = d->arr + j;
      s0 += ds->s; /* phi' = 1 */
      s1 += ds->sf * phi;
    }
    den = 2.0 * m * tot * dx;
    d->his[i] = s0/den;
    d->rho[i] = (s0 + s1)/den;
    d->jl[i] = i - j0;
    d->jr[i] = j1 - i;
  }
}

/* make fixed window */
INLINE void distr_winfix(distr_t *d, int m)
{
  int i, n = d->n;

  printf("fixed window, half width m = %d (%g)\n", m, m * d->dx);
  for (i = 0; i <= n; i++) {
    d->jl[i] = intmax(i - m, 0);
    d->jr[i] = intmin(i + m, n);
  }
}

/* make fixed window, guess width from 2 m dx = 1/sig(mf) */
INLINE void distr_winfixinvmf(distr_t *d)
{
  int i, m, n = d->n;
  double his, sig;

  /* error of mean force as variance divided by the number of samples */
  for (sig = his = 0, i = 0; i < n; i++) {
    distrsum_t *ds = d->arr + i;
    if (ds->s > 0.) {
      sig += ds->sf2 - ds->sf * ds->sf / ds->s;
      his += ds->s;
    }
  }
  sig = sqrt(sig/his);
  printf("mean force deviation is %g\n", sig);

  m = intmax(1, (int)( 0.5/(sig * d->dx) + .5 ));
  distr_winfix(d, m);
}

/* make adaptive window */
INLINE void distr_winadp(distr_t *d)
{
  int i, n = d->n;
  double his, sig, dx = d->dx, *hsum, *hhis;

  xnew(hsum, n + 1);
  xnew(hhis, n + 1);
  for (i = 0; i < n; i++)
    hhis[i] = .5 * (d->rho[i] + d->rho[i + 1]);
  for (hsum[0] = 0., i = 0; i < n; i++)
    hsum[i + 1] = hsum[i] + hhis[i];

  for (sig = his = 0, i = 0; i < n; i++) {
    distrsum_t *ds = d->arr + i;
    if (ds->s > 0.) {
      sig += ds->sf2 - ds->sf * ds->sf / ds->s;
      his += ds->s;
    }
  }
  sig = sqrt(sig/his);

  for (i = 0; i <= n; i++) {
    int k = 0, jl = i, jr = i+1, el = 0, er = 0;

    for (k = 0; el == 0 || er == 0; k++) {
      if (k % 2 == 0) { /* do the left */
        if (el || jl <= 1) {
          el = 1;
          continue;
        }
        if (hhis[jl - 1] / (hsum[jr] - hsum[jl] + 1e-10) < sig * dx) {
          el = 1;
        } else {
          jl--;
        }
      } else { /* do the right */
        if (er || jr >= n) {
          er = 1;
          continue;
        }
        if (hhis[jr] / (hsum[jr] - hsum[jl] + 1e-10) < sig * dx) {
          er = 1;
        } else {
          jr++;
        }
      }
    }
    d->jl[i] = jl;
    d->jr[i] = jr;
  }
  free(hhis);
  free(hsum);
}

/* obtain distribution from integrating the mean force
 * limiting case of distr_ii0 with infinite window size */
INLINE void distr_iimf(distr_t *d)
{
  int i, j, jl, jr, n = d->n;
  double x, max, den, tot, dx = d->dx;

  /* construct ln(rho) from the mean force */
  for (max = d->lnrho[0] = 0, i = 0; i < n; i++) {
    d->lnrho[i + 1] = d->lnrho[i] + d->mf[i] * dx;
    if (d->lnrho[i + 1] > max) max = d->lnrho[i + 1];
  }

  for (tot = 0, i = 0; i <= n; i++)
    tot += d->rho[i] = exp(d->lnrho[i] - max);
  tot *= dx;

  for (i = 0; i <= n; i++)
    d->rho[i] /= tot;
}

/* modulated integral identity from a fixed window [i - m, i + m),
 * set rho[], lnrho[], his[]
 * must call mfX(), first */
INLINE void distr_ii0(distr_t *d, int m)
{
  int i, j, jl, jr, n = d->n;
  double x, den, tot, dx = d->dx;

  /* construct ln(rho) from the mean force */
  for (i = 0; i < n; i++)
    d->lnrho[i + 1] = d->lnrho[i] + d->mf[i] * dx;

  for (d->hsum[0] = 0., i = 0; i < n; i++)
    d->hsum[i + 1] = d->hsum[i] + d->arr[i].s;
  tot = d->hsum[n]; /* total number of visits */
  printf("tot %g\n", tot);

  if (m > 0) distr_winfix(d, m);
  else if (m == 0) distr_winfixinvmf(d);
  else if (m == -1) distr_winadp(d);

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i */
    jl = d->jl[i];
    jr = d->jr[i];
    d->his[i] = d->hsum[jr] - d->hsum[jl];

    /* integrate over exp(lnrho[j] - lnrho[i])
     * Note: we can loop from j0 to j1, with little difference */
    for (den = 0., j = jl; j <= jr; j++) {
      x = exp(d->lnrho[j] - d->lnrho[i]);
      den += (j == jl || j == jr) ? x * .5 : x;
    }
    den *= tot * dx;
    d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
    d->his[i] /= tot * dx * (jr - jl);
  }
}

INLINE void distr_iiadp(distr_t *d)
{
  distr_iimf(d);
  distr_ii0(d, -1); /* use variable window */
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
  xnew(d->his, d->n + 1);
  xnew(d->mf, d->n + 1);
  xnew(d->mdv, d->n + 1);
  xnew(d->jl, d->n + 1);
  xnew(d->jr, d->n + 1);
  xnew(d->err, d->n + 1);
  xnew(d->hsum, d->n + 1);
  return d;
}

INLINE void distr_close(distr_t *d)
{
  if (!d) return;
  free(d->arr);
  free(d->rho);
  free(d->lnrho);
  free(d->his);
  free(d->mf);
  free(d->mdv);
  free(d->jl);
  free(d->jr);
  free(d->err);
  free(d->hsum);
  free(d);
}

INLINE void distr_add(distr_t *d, double x, double f, double dv, double w)
{
  if (x >= d->xmin && x <= d->xmax) {
    distrsum_t *ds = d->arr + (int)((x - d->xmin)/d->dx);
    ds->s += w;
    ds->sf += w * f;
    ds->sf2 += w * f * f;
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
  fprintf(fp, "# 1 %d %.14e %.14e\n", d->n, d->xmin, d->dx);

  for (i = i0; i <= i1; i++) {
    double f = 0.0, f2 = 0.0, dv = 0.0;
    ds = d->arr + i;
    if (i < n && ds->s > 0.) {
      f = ds->sf / ds->s;
      f2 = ds->sf2 / ds->s - f * f; /* convert to variance */
      dv = ds->sdv / ds->s;
    }
    fprintf(fp, "%g %.14e %.14e %.14e %.14e "
        "%.14e %.14e %.14e %.14e %.14e %g %g %g\n",
      d->xmin + i * d->dx, ds->s, f, f2, dv,
      d->rho[i], d->lnrho[i], d->his[i], d->mf[i], d->mdv[i],
      d->xmin + d->dx * d->jl[i], d->xmin + d->dx * d->jr[i], d->err[i]);
  }
  fclose(fp);
  return 0;
}

/* load a previous histogram */
STRCLS int distr_load(distr_t *d, const char *fn)
{
  FILE *fp;
  char s[1024];
  int ver, i, n = d->n, nlin;
  double x, sm, f, f2, dv, xmin = d->xmin, dx = d->dx;
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
    if (5 != sscanf(s, "%lf%lf%lf%lf%lf", &x, &sm, &f, &f2, &dv)) {
      fprintf(stderr, "sscanf error\n");
      goto ERR;
    }
    /* locate the bin */
    i = (int)( (x - xmin)/dx + .1);
    if (i < 0 || i >= n) {
      fprintf(stderr, "bad x %g, xmin %g, i %d, n %d\n", x, xmin, i, n);
      goto ERR;
    }
    ds = d->arr + i;
    ds->s = sm;
    ds->sf = f * sm;
    ds->sf2 = (f2 + f*f)*sm;
    ds->sdv = dv * sm;
    if (i >= n - 1) break;
  }
  return 0;
ERR:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  /* clear everything on error */
  for (i = 0; i <= n; i++) distrsum_clear(d->arr + i);
  return -1;
}

#endif
