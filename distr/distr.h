#include "def.h"
#include "util.h"
#include "specfunc.h"
#include "rc.h"
#ifndef DISTR_H__
#define DISTR_H__

#define DSDIF_HIS1 0x01 /* the first distribution is histogram */
#define DSDIF_HIS2 0x02 /* the second (reference) is histogram */

/* return KS difference between f and ref, with p-value
 * f and ref are distribution densities
 * if ishis & DSDIF_HIS1, f[0..n-1] is histogram, otherwise f[0..n] is the distribution
 * is ishis & DSDIF_HIS2, ref[0..n-1] is histogram, otherwise ref[0..n] is the distribution
 * smpsize is the sample size */
INLINE double ksdif1d(const double f[], const double ref[], int n, int ishis,
    double smpsiz, double *p, double err[])
{
  int i, alloc = 0, his1 = ishis & DSDIF_HIS1, his2 = ishis & DSDIF_HIS2;
  double tot1, tot2, sum1, sum2, e, emax;
 
  if (err == NULL) { alloc = 1; xnew(err, n + 1); }
  err[0] = 0.; 

  /* compute the normalization factor for both distributions */
  for (tot1 = tot2 = 0., i = 0; i < n; i++) {
    tot1 += his1 ? f[i]   : .5*(f[i] + f[i+1]);
    tot2 += his2 ? ref[i] : .5*(ref[i] + ref[i+1]);
  }
  if (tot1 <= 0. || tot2 <= 0.) {
    fprintf(stderr, "ksdif: bad distribution, tot 1: %g, 2: %g\n", tot1, tot2);
    return 0;
  }

  for (sum1 = sum2 = 0., emax = 0., i = 0; i < n; i++) {
    sum1 += (his1 ? f[i]   : .5*(f[i] + f[i+1])) / tot1;
    sum2 += (his2 ? ref[i] : .5*(ref[i] + ref[i+1])) / tot2;
    err[i + 1] = e = sum1 - sum2;
    if ((e = fabs(e)) > emax) emax = e;
  }
  e = sqrt(smpsiz);
  emax *= e + .12 + .11/e; /* a formula from numerical recipes */
  if (p != NULL) *p = ksq(emax);
  if (alloc) free(err);
  return emax;
}

/* entropic difference = integral log(ref/r) ref(x) dx */
INLINE double entdif1d(const double f[], const double ref[], int n,
    int ishis, double err[])
{
  int i, alloc = 0, his1 = ishis & DSDIF_HIS1, his2 = ishis & DSDIF_HIS2;
  double tot1, tot2, y1, y2, e;
 
  if (err == NULL) { alloc = 1; xnew(err, n + 1); }
  err[0] = 0.; 

  /* compute the normalization factor for both distributions */
  for (tot1 = tot2 = 0., i = 0; i < n; i++) {
    tot1 += his1 ? f[i]   : .5*(f[i] + f[i+1]);
    tot2 += his2 ? ref[i] : .5*(ref[i] + ref[i+1]);
  }
  if (tot1 <= 0. || tot2 <= 0.) {
    fprintf(stderr, "entdif: bad distribution, tot 1: %g, 2: %g\n", tot1, tot2);
    return 0;
  }

  for (e = 0., i = 0; i < n; i++) {
    y1 = (his1 ? f[i]   : .5*(f[i] + f[i+1])) / tot1;
    y2 = (his2 ? ref[i] : .5*(ref[i] + ref[i+1])) / tot2;
    if (y1 <= 0 || y2 <= 0) {
      if (err) err[i] = 0;
      continue;
    }
    y2 = log(y1/y2);
    if (err) err[i] = y2;
    e += y2 * y1;
  }
  if (alloc) free(err);
  return e;
}


typedef struct {
  double s;
  double sf;  /* mean force */
  double sf2; /* mean force^2 */
  double sdv; /* average value for the 2nd derivative */
} distrsum_t;

INLINE void distrsum_clear(distrsum_t *x) { x->s = x->sf = x->sf2 = x->sdv = 0.0; }

typedef struct { int il, ir; } distrwin_t;

typedef struct {
  int n;
  double xmin, xmax, dx;
  distrsum_t *arr;
  double *rho, *lnrho, *his;
  double *mf, *mdv;
  distrwin_t *win; /* window size */
  double *err;
  double *hsum; /* histogram sum for comparison */
  double *mfsig; /* std. dev. of mf, for adaptive window */
  double *norm; /* normalization vector */
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
 * linear phi(x) is used for the integral identity
 * improves mean force, slightly smoothness of distribution itself */
INLINE void distr_mfwin(distr_t *d, int m, int ii)
{
  int i, j, j0, j1, n = d->n;
  double s0, s1, s, phi, dx = d->dx;

  distr_mf0(d);
  for (i = 0; i < n; i++) {
    j0 = intmax(0, i - m);
    j1 = intmin(n, i + m + 1);
    for (s0 = s1 = s = 0., j = j0; j < j1; j++) {
      phi = (j == i) ? 0 : (j + .5 - (j < i ? j0 : j1)) * dx;
      s0 += d->mf[j];
      s1 += d->mdv[j] * phi;
      s += 1;
    }
    if (ii) s0 += s1; /* add correction from the second-order derivative */
    d->mf[i] = (s > 0) ? s0/s : 0.0;
  }
}

/* Adib-Jazynski identity from [i - m, i + m)
 * linear function phi(x) is assumed */
INLINE void distr_aj(distr_t *d, int m)
{
  int i, j, j0, j1, n = d->n;
  double s0, s1, phi, den, tot, s, sf, dx = d->dx;
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
      s = ds->s;
      sf = ds->sf;
      if (d->norm) {
        double nm = .5 * (d->norm[j] + d->norm[j + 1]);
        s /= nm;
        sf /= nm;
      }
      s0 += s; /* phi' = 1 */
      s1 += sf * phi;
    }
    den = 2.0 * m * tot * dx;
    d->his[i] = s0/den;
    d->rho[i] = (s0 + s1 > 0) ? (s0 + s1)/den : 0.;
    d->win[i].il = i - j0;
    d->win[i].ir = j1 - i;
  }
}

/* make fixed windows
 * if m == 0, guess width from 2 m dx = 1/sig(mf) */
INLINE void distr_winfixed(distr_t *d, distrwin_t *win, double gam, int m)
{
  int i, n = d->n;
  double his, sig;

  if (m == 0) { /* guess window size */
    /* error of mean force as variance divided by the number of samples */
    for (sig = his = 0, i = 0; i < n; i++) {
      distrsum_t *ds = d->arr + i;
      if (ds->s > 0.) {
        sig += ds->sf2 - ds->sf * ds->sf / ds->s;
        his += ds->s;
      }
    }
    sig = sqrt(sig/his);

    m = intmax(1, (int) (0.5 * gam/(sig * d->dx) + .5));
    printf("mean force standard deviation is %g, m %d\n", sig, m);
  }

  for (i = 0; i <= n; i++) {
    win[i].il = intmax(i - m, 0);
    win[i].ir = intmin(i + m, n);
  }  
}

/* estimate local standard deviation of the mean force
 * average over local window of 2 m */
INLINE void distr_mfsig0(distr_t *d, double sampmin)
{
  int i, m, n = d->n;
  double x, sig, his, *asig, *hsum;

  xnew(asig, n + 1);
  xnew(hsum, n + 1);
  for (asig[0] = hsum[0] = 0, i = 0; i < n; i++) {
    distrsum_t *ds = d->arr + i;
    hsum[i+1] = hsum[i] + ds->s;
    x = (ds->s > 0) ? (ds->sf2 - ds->sf * ds->sf / ds->s) : 0.;
    asig[i+1] = asig[i] + x;
  }
  for (i = 0; i < n; i++) {
    /* start from the minimal window size, expand until included sampmin data points */
    for (m = 1; m < n; m++) {
      int j0 = intmax(0, i - m), j1 = intmin(i + m + 1, n);
      his = hsum[j1] - hsum[j0];
      sig = asig[j1] - asig[j0];
      d->mfsig[i] = his > 0. ? sqrt(sig/his) : 0.;
      if (his > sampmin || (j0 == 0 && j1 == n)) break;
    }
  }
  d->mfsig[n] = d->mfsig[n-1];
  free(asig);
  free(hsum);
} 

/* make symmetrical but x-dependent window
 * requires d->mfsig */
INLINE void distr_winadp0(distr_t *d, distrwin_t *win, double gam, int mlimit)
{
  int i, m, n = d->n;
  double dx = d->dx;

  if (mlimit < 0) mlimit = 100000000; /* set a very large number */

  for (i = 0; i <= n; i++) {
    m = (int) intmax(1, (int) (0.5 * gam/(d->mfsig[i] * dx) + .5));
    win[i].il = intmax(0, i - m);
    win[i].ir = intmin(n, i + m);
  }
}

/* make adaptive window
 * requires d->rho and d->mfsig */
INLINE void distr_winadp(distr_t *d, int mlimit)
{
  int i, n = d->n;
  double dx = d->dx, *hsum, *hhis;

  xnew(hsum, n + 1);
  xnew(hhis, n + 1);
  for (i = 0; i < n; i++) {
    hhis[i] = .5 * (d->rho[i] + d->rho[i + 1]);
    if (d->norm) hhis[i] *= .5 * (d->norm[i] + d->norm[i + 1]);
  }
  for (hsum[0] = 0., i = 0; i < n; i++)
    hsum[i + 1] = hsum[i] + hhis[i];

  if (mlimit < 0) mlimit = 100000000; /* set a very large number */

  for (i = 0; i <= n; i++) {
    int k = 0, jl = i, jr = i+1, el = 0, er = 0, jlm, jrm;

    jlm = intmax(0, i - mlimit);
    jrm = intmin(n, i + mlimit);
    for (k = 0; el == 0 || er == 0; k++) {
      if (k % 2 == 0) { /* do the left */
        if (el || jl <= jlm) {
          el = 1;
          continue;
        }
        if (hhis[jl - 1] / (hsum[jr] - hsum[jl] + 1e-10) < d->mfsig[jl - 1] * dx) {
          el = 1;
        } else {
          jl--;
        }
      } else { /* do the right */
        if (er || jr >= jrm) {
          er = 1;
          continue;
        }
        if (hhis[jr] / (hsum[jr] - hsum[jl] + 1e-10) < d->mfsig[jr] * dx) {
          er = 1;
        } else {
          jr++;
        }
      }
    }
    d->win[i].il = jl;
    d->win[i].ir = jr;
  }
  free(hhis);
  free(hsum);
}

/* obtain distribution from integrating the mean force
 * limiting case of distr_ii0 with infinite window size */
INLINE void distr_iimf(distr_t *d)
{
  int i, n = d->n;
  double max, tot, x, lntot, dx = d->dx;

  /* integrate mf to get ln(rho) */
  for (max = d->lnrho[0] = 0, i = 0; i < n; i++) {
    d->lnrho[i + 1] = d->lnrho[i] + d->mf[i] * dx;
    if (d->lnrho[i + 1] > max) max = d->lnrho[i + 1];
  }

  for (tot = 0, i = 0; i <= n; i++) {
    x = exp(d->lnrho[i] -= max);
    if (d->norm) x *= d->norm[i];
    tot += d->rho[i] = x;
  }
  tot *= dx;
  lntot = log(tot);

  for (i = 0; i <= n; i++) {
    d->lnrho[i] -= lntot;
    d->rho[i] /= tot;
  }
}

/* modulated integral identity
 * set rho[], lnrho[], his[]
 * before calling this function
 * first call mfX() to compute the mean force, 
 * and call winX() to compute the window */
INLINE void distr_ii0(distr_t *d, distrwin_t *win)
{
  int i, j, jl, jr, n = d->n;
  double x, den, tot, dx = d->dx;

  /* construct ln(rho) from the mean force */
  distr_iimf(d);

  for (d->hsum[0] = 0., i = 0; i < n; i++)
    d->hsum[i + 1] = d->hsum[i] + d->arr[i].s;
  tot = d->hsum[n]; /* total number of visits */
  
  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i */
    jl = win[i].il;
    jr = win[i].ir;
    d->his[i] = d->hsum[jr] - d->hsum[jl];

    /* integrate over exp(lnrho[j] - lnrho[i])
     * Note: we can loop from j0 to j1, with little difference */
    for (den = 0., j = jl; j <= jr; j++) {
      x = exp(d->lnrho[j] - d->lnrho[i]);
      if (d->norm) x *= d->norm[j];
      den += (j == jl || j == jr) ? x * .5 : x;
    }
    den *= tot * dx;
    d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
    d->his[i] /= tot * dx * (jr - jl);
  }
}

/* integral identity with an adaptive window
 * may not always help, use with care */
INLINE void distr_iiadp(distr_t *d, int mlimit, double sigsampmin)
{
  /* construct a tentative d->rho from integrating the mean force */
  distr_iimf(d);

  /* estimate mean force variance */
  distr_mfsig0(d, sigsampmin);

  /* compute adaptive windows */
  distr_winadp(d, mlimit);

  /* estimate the distribution */
  distr_ii0(d, d->win);
}

/* perform integral identity
 * iitype: 0: Adib-Jazynski, 1: fraction identity, 2: integrate the mean force
 * halfwin: half window size. > 0: fixed size, 0: fixed size (determined
 *          automatically), -1: adaptive window size
 * mfhalfwin: half window size for mean force integration
 * gamma: window size amplification factor
 * mlimit: maximal half window size (-1 to disable)
 * sampmin: minimal number of samples in estimating mean force std. */
INLINE void distr_iiez(distr_t *d, int iitype, int halfwin, int mfhalfwin,
    double gam, int mlimit, double sampmin)
{
  if (iitype == 0) {
    distr_aj(d, halfwin);
  } else if (iitype == 1 || iitype == 2) {
    /* compute mean force */
    if (mfhalfwin == 0) {
      distr_mf0(d);
    } else {
      distr_mfwin(d, abs(mfhalfwin), mfhalfwin > 0 ? 1 : 0);
    }

    /* compute the distribution */
    if (iitype == 1) {
      if (halfwin >= -1) {
        if (halfwin >= 0) {
          distr_winfixed(d, d->win, gam, halfwin);
        } else {
          distr_mfsig0(d, sampmin);
          distr_winadp0(d, d->win, gam, mlimit);
        }
        distr_ii0(d, d->win);
      } else {
        distr_iiadp(d, mlimit, sampmin);
      }
    } else {
      distr_iimf(d);
    }
  }
}

INLINE distr_t *distr_open(double xmin, double xmax, double dx)
{
  distr_t *d;
  xnew(d, 1);
  d->xmin = xmin;
  d->dx = dx;
  d->n = (int)((xmax - xmin)/dx + 0.999999f);
  d->xmax = xmin + d->n * dx;
  xnew(d->arr, d->n + 1);
  xnew(d->rho, d->n + 1);
  xnew(d->lnrho, d->n + 1);
  xnew(d->his, d->n + 1);
  xnew(d->mf, d->n + 1);
  xnew(d->mdv, d->n + 1);
  xnew(d->win, d->n + 1);
  xnew(d->err, d->n + 1);
  xnew(d->hsum, d->n + 1);
  xnew(d->mfsig, d->n + 1);
  d->norm = NULL;
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
  free(d->win);
  free(d->err);
  free(d->hsum);
  free(d->mfsig);
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
        "%.14e %.14e %.14e %.14e %.14e %g %g "
        "%g %g\n",
      d->xmin + i * d->dx, ds->s, f, f2, dv,
      d->rho[i], d->lnrho[i], d->his[i], d->mf[i], d->mdv[i],
      d->xmin + d->dx * d->win[i].il, d->xmin + d->dx * d->win[i].ir, 
      d->err[i], d->mfsig[i]);
  }
  fclose(fp);
  return 0;
}

/* load a previous histogram */
INLINE int distr_load(distr_t *d, const char *fn)
{
  FILE *fp;
  char s[1024];
  int ver, i, n = d->n, nlin;
  double x, sm, f, f2, dv, rho, xmin = d->xmin, dx = d->dx;
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
    if (6 != sscanf(s, "%lf%lf%lf%lf%lf%lf", &x, &sm, &f, &f2, &dv, &rho)) {
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
    d->rho[i] = rho;
    if (i >= n - 1) break;
  }
  return 0;
ERR:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  /* clear everything on error */
  for (i = 0; i <= n; i++) distrsum_clear(d->arr + i);
  return -1;
}

/* least square solver for a 2D potential from the mean forces
 * assume periodic along both directions */
typedef struct {
  int n, m; /* dimension of x and y */
  real dx, dy;
  rcomplex_t *u, *f, *g; /* pmf, mean force along x and y */
  rcomplex_t *ru, *rf, *rg, *tmp; /* reciprocal space  */
  rcomplex_t *en, *em; /* exp(i 2 pi I / n) and exp(j 2 pi I / m) */
} ftsolver2d_t;

INLINE ftsolver2d_t *ftsolver2d_open(int n, real dx, int m, real dy)
{
  ftsolver2d_t *fts;
  int i;

  xnew(fts, 1);
  fts->n = n;
  fts->m = m;
  fts->dx = dx;
  fts->dy = dy;
  xnew(fts->u, n * m);
  xnew(fts->f, n * m);
  xnew(fts->g, n * m);
  xnew(fts->ru, n * m);
  xnew(fts->rf, n * m);
  xnew(fts->rg, n * m);
  xnew(fts->tmp, n * m);
  xnew(fts->en, 2 * n);
  for (i = 0; i < 2 * n; i++) {
    double phi = M_PI*i/n;
    fts->en[i].re = (real) cos(phi);
    fts->en[i].im = (real) sin(phi);
  }
  xnew(fts->em, 2 * m);
  for (i = 0; i < 2 * m; i ++) {
    double phi = M_PI*i/m;
    fts->em[i].re = (real) cos(phi);
    fts->em[i].im = (real) sin(phi);
  }
  return fts;
}

INLINE void ftsolver2d_close(ftsolver2d_t *fts)
{
  free(fts->u);
  free(fts->f);
  free(fts->g);
  free(fts->ru);
  free(fts->rf);
  free(fts->rg);
  free(fts->tmp);
  free(fts->en);
  free(fts->em);
  free(fts);
}

/* slow Fourier transform rf of f
 * decomposed along x and y 
 * sgn:   0 to get coeffients, 1 to combine components
 * inpre: 1 if input is real */
INLINE int ftsolver2d_ft(ftsolver2d_t *fts, rcomplex_t *rf, const rcomplex_t *f,
   unsigned flags)
{
  int i, j, k, l, n = fts->n, m = fts->m;
  int sgn = (flags & 0x1), inpre = (flags & 0x2);
  rcomplex_t fkl, x;

  /* FT along the m direction, save to tmp */
  for (i = 0; i < n; i++) { /* for different rows */
    for (l = 0; l < m; l++) { /* for different wave vectors */
      for (fkl = rc_zero, j = 0; j < m; j++) { /* integrate */
        x = rc_mul(fts->em[(j * l) % m * 2], f[i * m + j]);
        fkl = rc_add(fkl, x);
      }
      fts->tmp[i*m + l] = fkl;
      if (inpre) { /* reflect around l = m/2 */
        if (l >= m/2) break;
        if (l > 0 && l*2 < m) fts->tmp[i*m + m - l] = rc_star(fkl);
      }
    }
  }

  /* FT along the n direction, save to rf */
  for (l = 0; l < m; l++) { /* for different columns */
    for (k = 0; k < n; k++) { /* for different wave vectors */
      for (fkl = rc_zero, i = 0; i < n; i++) { /* integrate */
        x = rc_mul(fts->en[(i * k) % n * 2], fts->tmp[i * m + l]);
        fkl = rc_add(fkl, x);
      }
      if (!sgn) fkl = rc_sdiv(rc_star(fkl), (real)(n * m)); /* fkl* / (n m) */
      rf[k*m + l] = fkl;
      if (inpre && l > 0 && l*2 < m) /* reflect: rf(n - k, m - l) = rf(n, l)* */
        rf[(n - k) % n * m + (m - l)] = rc_star(fkl);
    }
    if (inpre && l >= m/2) break;
  }
  return 0;
}

/* return the least-square 2D ln(rho)
 * from the mean forces along the two directions 
 * fx and fy are (n+1)*(m+1) */
INLINE void ftsolver2d_solve(ftsolver2d_t *fts, double *u, double *fx, double *fy)
{
  int i, j, id, id1, k, l, kl, n = fts->n, m = fts->m;
  real ck, sk, cl, sl, c2k, c2l, dx = fts->dx, dy = fts->dy;
  rcomplex_t x;

  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    id = i * m + j;
    id1 = i * (m + 1) + j;
    fts->f[id].re = (real) fx[id1];
    fts->f[id].im = 0.f;
    fts->g[id].re = (real) fy[id1];
    fts->g[id].im = 0.f;
  }
  ftsolver2d_ft(fts, fts->rf, fts->f, 0x2);
  ftsolver2d_ft(fts, fts->rg, fts->g, 0x2);
  /* the loop is O(nm), no need to optimize */
  for (k = 0; k < n; k++) {
    ck = fts->en[k].re;
    sk = fts->en[k].im;
    c2k = fts->en[k*2].re;
    for (l = 0; l < m; l++) {
      kl = k*m + l;
      if (kl == 0) {
        fts->ru[kl].re = fts->ru[kl].im = 0.f;
        continue;
      }
      cl = fts->em[l].re;
      sl = fts->em[l].im;
      c2l = fts->em[l*2].re;
      /* sin(pi k/N) cos(pi l/M) f(k, l) dx */
      x = rc_smul(fts->rf[kl], sk * cl * dx);
      /* + cos(pi k/N) sin(pi l/M) g(k, l) dy */
      x = rc_add(x, rc_smul(fts->rg[kl], ck * sl * dy));
      /* multiply -i exp[ - pi i (k/N + l/M) ] */
      x = rc_mul(x, rc_make(-(sk*cl + ck*sl), -(ck*cl - sk*sl)));
      if (fabs(c2k*c2l - 1) < 1e-15) {
        x.re = x.im = 0.f;
      } else {
        x = rc_sdiv(x, 1 - c2k*c2l);
      }
      fts->ru[kl] = x;
    }
  }
  ftsolver2d_ft(fts, fts->u, fts->ru, 1); /* inverse transform back */
  if (u != NULL) {
    for (i = 0; i <= n; i++)
      for (j = 0; j <= m; j++) {
        id = (i%n) * m + (j%m);
        id1 = i * (m + 1) + j;
        u[id1] = fts->u[id].re;
      }
  }
}

typedef struct { double s, sf, sg, sf2, sg2, sfg; } distr2dsum_t;

INLINE void distr2dsum_clear(distr2dsum_t *x)
{ x->s = x->sf = x->sg = x->sf2 = x->sg2 = 0.0; }

typedef struct { int il, ir, jl, jr; } distr2dwin_t;

typedef struct {
  int n, m;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
  distr2dsum_t *arr;
  double *rho, *lnrho, *his;
  double *mf, *mg;
  double tot;
  distr2dwin_t *win;
} distr2d_t;

INLINE distr2d_t *distr2d_open(double xmin, double xmax, double dx,
    double ymin, double ymax, double dy)
{
  distr2d_t *d;
  int n1m1;

  xnew(d, 1);
  
  d->xmin = xmin;
  d->dx = dx;
  d->n = (int)((xmax - xmin)/dx + 0.5f);
  d->xmax = xmin + d->n * dx;

  d->ymin = ymin;
  d->dy = dy;
  d->m = (int)((ymax - ymin)/dy + 0.5f);
  d->ymax = ymin + d->m * dy;
  d->tot = 0;

  n1m1 = (d->n + 1) * (d->m + 1);
  xnew(d->arr, n1m1);
  xnew(d->rho, n1m1);
  xnew(d->lnrho, n1m1);
  xnew(d->his, n1m1);
  xnew(d->mf, n1m1);
  xnew(d->mg, n1m1);
  xnew(d->win, n1m1);
  return d;
}

INLINE void distr2d_close(distr2d_t *d)
{
  if (!d) return;
  free(d->arr);
  free(d->rho);
  free(d->lnrho);
  free(d->his);
  free(d->mf);
  free(d->mg);
  free(d->win);
  free(d);
}

INLINE void distr2d_add(distr2d_t *d, double x, double y, double f, double g, double w)
{
  if (x >= d->xmin && x <= d->xmax && y >= d->ymin && y <= d->ymax) {
    int i = (int)((x - d->xmin)/d->dx), j = (int)((y - d->ymin)/d->dy);
    int m1 = d->m + 1;
    distr2dsum_t *ds = d->arr + i * m1 + j; 
    ds->s   += w;
    ds->sf  += w * f;
    ds->sg  += w * g;
    ds->sf2 += w * f * f;
    ds->sg2 += w * g * g;
    ds->sfg += w * f * g;
  }
}

INLINE double distr2d_tot(distr2d_t *d)
{
  int i, j, n = d->n, m = d->m;
  double tot = 0;
  
  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    tot += d->arr[i * (m + 1) + j].s;
  }
  return d->tot = tot;
}

/* compute the mean force from a single bin */
INLINE void distr2d_mf0(distr2d_t *d)
{
  int i, j, id, n = d->n, m = d->m, m1 = d->m + 1;
  distr2dsum_t *ds;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      double mf = 0., mg = 0.;
      id = i * m1 + j;
      ds = d->arr + id;
      if (ds->s > 0.) {
        mf = ds->sf / ds->s;
        mg = ds->sg / ds->s;
      }
      d->mf[id] = mf;
      d->mg[id] = mg;
    }
}

/* compute the mean force from a window at least of size 2w+1 */
INLINE void distr2d_mfwin(distr2d_t *d, int w)
{
  int i, j, id, n = d->n, m = d->m, m1 = d->m + 1, ww;
  distr2dsum_t *ds;

  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    double ms, mf, mg;
    /* try to expand the window until we have something */
    for (ww = w; ww < intmin(n/2, m/2); ww++) {
      int il = i - ww, ir = i + ww + 1, jl = j - ww, jr = j + ww + 1;
      int i1, j1;
      ms = mf = mg = 0.;
      for (i1 = il; i1 < ir; i1++)
      for (j1 = jl; j1 < jr; j1++) {
        id = (i1 + n) % n * m1 + (j1 + m) % m;
        ds = d->arr + id;
        if (ds->s <= 0) continue;
        ms += ds->s;
        mf += ds->sf;
        mg += ds->sg;
      }
      if (ms > 0) break; 
    }
    id = i * m1 + j;
    d->mf[id] = (ms > 0) ? mf/ms : 0;
    d->mg[id] = (ms > 0) ? mg/ms : 0;
  }
}

/* compute lnrho from least-square fitting of the mean force */
INLINE void distr2d_calclnrho(distr2d_t *d)
{
  int i, j, n = d->n, m = d->m;
  double tot = 0, lntot;
  
  ftsolver2d_t *fts = ftsolver2d_open(d->n, (real) d->dx, d->m, (real) d->dy);
  ftsolver2d_solve(fts, d->lnrho, d->mf, d->mg);
  ftsolver2d_close(fts);

  /* normalize d->lnrho */
  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    tot += exp(d->lnrho[i * (m + 1) + j]);
  }
  tot *= d->dx * d->dy;
  lntot = log(tot);
  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++)
    d->lnrho[i * (m + 1) + j] -= lntot;
}

/* compute the mean force variance */
INLINE double distr2d_mfvar(distr2d_t *d)
{
  int i, j, n = d->n, m = d->m;
  double sm = 0, f2sm = 0, s, f, g;

  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    distr2dsum_t *ds = d->arr + i * (m + 1) + j;
    if ((s = ds->s) <= 0.) continue;
    sm += 1;
    f = ds->sf / s;
    g = ds->sg / s;
    f2sm += .5*(ds->sf2/s - f*f + ds->sg2/s - g*g);
  }
  return sm > 0. ? f2sm/sm : 0.;
}

/* fixed window */
INLINE void distr2d_winfixed(distr2d_t *d, distr2dwin_t *dwin, int w, double smin)
{
  int i, j, n = d->n, m = d->m, m1 = d->m + 1, ww, ww1;
  distr2dsum_t *ds;

  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    distr2dwin_t *win = dwin + i * m1 + j;
    win->il = i - w;
    win->ir = i + w;
    win->jl = j - w;
    win->jr = j + w;

    if (smin > 0) {
      for (ww = w; ww < intmin(n/2, m/2); ww++) {
        int il = i - ww, ir = i + ww, jl = j - ww, jr = j + ww, i1, j1;
        double ms = 0.;
        for (i1 = il; i1 < ir; i1++)
        for (j1 = jl; j1 < jr; j1++) {
          ms += d->arr [ (i1 + n) % n * m1 + (j1 + m) % m ].s;
        }
        if (ms > 0) break; 
      }
      win->il = i - ww;
      win->ir = i + ww;
      win->jl = j - ww;
      win->jr = j + ww;
    }
  }
}

/* integral identity, assume periodic function */
INLINE void distr2d_ii0(distr2d_t *d, distr2dwin_t *win)
{
  int i, j, id, il, ir, jl, jr, ii, jj, iip, jjp, ij, n = d->n, m = d->m;
  double num, den, den0, wt, dlnrho, tot;

  tot = distr2d_tot(d);
  for (i = 0; i < n; i++)
  for (j = 0; j < m; j++) {
    id = i * (m + 1) + j;
    il = win[id].il;
    ir = win[id].ir;
    jl = win[id].jl;
    jr = win[id].jr;

    num = den = den0 = 0.;
    /* get the number of visits to the region */
    for (ii = il; ii < ir; ii++)
    for (jj = jl; jj < jr; jj++) {
      iip = (ii + n) % n;
      jjp = (jj + n) % n;
      num += d->arr[iip * (m + 1) + jjp].s;
    }
    for (ii = il; ii <= ir; ii++)
    for (jj = jl; jj <= jr; jj++) {
      iip = (ii + n) % n;
      jjp = (jj + n) % n;
      ij = iip * (m + 1) + jjp;
      wt = 1.0;
      if (jj == jl || jj == jr) wt *= .5;
      if (ii == il || ii == ir) wt *= .5;
      den0 += wt;
      dlnrho = d->lnrho[ij] - d->lnrho[id];
      den += wt * exp(dlnrho);
    }
    den *= tot * d->dx * d->dy;
    den0 *= tot * d->dx * d->dy;
    d->rho[id] = num/den; /* unbiased result */
    d->his[id] = num/den0; /* biased result */
  }
}

INLINE int distr2d_save(distr2d_t *d, const char *fn)
{
  FILE *fp;
  int i, i0 = 0, i1 = d->n, n = d->n, id;
  int j, j0 = 0, j1 = d->m, m = d->m, m1 = d->m + 1;
  distr2dsum_t *ds = d->arr;

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# 0 %d %.14e %.14e %d %.14e %.14e\n",
      d->n, d->xmin, d->dx, d->m, d->ymin, d->dy);

  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {
      double sm = 0., f = 0.0, g = 0.0, f2 = 0.0, g2 = 0.0, fg = 0.0;
      id = i * m1 + j;
      if (i < n && j < m) {
        ds = d->arr + id;
        sm = ds->s;
        if (i < n && j < m && sm > 0.) {
          f = ds->sf / ds->s;
          g = ds->sg / ds->s;
          f2 = ds->sf2 / ds->s - f * f; /* convert to variance */
          g2 = ds->sg2 / ds->s - f * f; /* convert to variance */
          fg = ds->sfg / ds->s - f * g; /* covariance */
        }
      }
      fprintf(fp, "%g %g %.14e %.14e %.14e %.14e %.14e %.14e %g %g %g\n",
        d->xmin + i * d->dx, d->ymin + j * d->dy, 
        sm, f, g, f2, g2, fg,
        d->rho[id], d->lnrho[id], d->his[id]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* load a previous histogram */
INLINE int distr2d_load(distr2d_t *d, const char *fn)
{
  FILE *fp;
  char s[1024];
  int ver, i, j, n = d->n, m = d->m, nlin;
  double x, y, sm, f, g, f2, g2, fg;
  double xmin = d->xmin, dx = d->dx, ymin = d->ymin, dy = d->dy;
  distr2dsum_t *ds;

  xfopen(fp, fn, "r", return -1);

  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  if (7 != sscanf(s, " # %d%d%lf%lf%d%lf%lf", &ver, &i, &x, &f, &j, &y, &g)
      || i != n || fabs(x - xmin) > 1e-3 || fabs(f - dx) > 1e-5
      || j != m || fabs(y - ymin) > 1e-3 || fabs(g - dy) > 1e-5) {
    fprintf(stderr, "error: n %d, %d; xmin %g, %g; dx %g, %g; "
        "m %d, %d; ymin %g, %g; dy %g, %g\n",
        i, n, x, xmin, f, dx, j, m, y, ymin, g, dy);
    fclose(fp);
    return -1;
  }

  /* only read in raw data */
  for (nlin = 2; fgets(s, sizeof s, fp); nlin++) {
    strip(s);
    if (s[0] == '\0') continue;
    if (8 != sscanf(s, "%lf%lf%lf%lf%lf%lf%lf%lf", 
          &x, &y, &sm, &f, &g, &f2, &g2, &fg)) { /* only scan raw data */
      fprintf(stderr, "sscanf error\n");
      goto ERR;
    }
    /* locate the bin */
    i = (int)( (x - xmin)/dx + .1);
    j = (int)( (y - ymin)/dy + .1);
    if (i < 0 || i > n || j < 0 || j > m) {
      fprintf(stderr, "bad x %g, xmin %g, i %d/%d, y %g, ymin %g, j %d/%d\n", 
          x, xmin, i, n, y, ymin, j, m);
      goto ERR;
    }
    if (i == n || j == m) continue;
    ds = d->arr + i * (m + 1) + j;
    ds->s = sm;
    ds->sf = f * sm;
    ds->sg = g * sm;
    ds->sf2 = (f2 + f*f)*sm;
    ds->sg2 = (g2 + g*g)*sm;
    ds->sfg = (fg + f*g)*sm;
  }
  return 0;
ERR:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  /* clear everything on error */
  for (i = 0; i <= (n + 1)*(m + 1); i++) distr2dsum_clear(d->arr + i);
  return -1;
}

#endif
