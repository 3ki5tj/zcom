#include "def.h"
#include "util.h"
#include "specfunc.c"
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



/* data collected in each bin */
typedef struct {
  double s;
  double sf;  /* mean force */
  double sf2; /* mean force^2 */
  double sdv; /* average value for the 2nd derivative */
} distrsum_t;

INLINE void distrsum_clear(distrsum_t *x) { x->s = x->sf = x->sf2 = x->sdv = 0.0; }


/* 1D window boundaries */
typedef struct { int il, ir; } distrwin_t;


/* distribution object */
typedef struct {
  int n;
  double xmin, xmax, dx;
  distrsum_t *arr; /* [0..n), bin data */
  double *rho; /* [0..n], output distribution */
  double *lnrho;  /* [0..n], from mean force integration */
  double *his; /* [0..n], histogram average */
  double *mf, *mdv; /* estimated mean force and deterministic part of mdv 
                     * in the simplest case, mf[i] = arr[i].sf/arr[i].s */
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
 * where `i' and `m' are integer bin index
 * linear function phi(x) is assumed
 * asymmetric windows are used at the boundaries
 * uses d->arr[j].s and .sf (but not the second derivatives .sf2) 
 * sets d->his, d->rho, d->win */
INLINE void distr_aj(distr_t *d, int m)
{
  int i, j, j0, j1, n = d->n;
  double s0, s1, phi, den, tot, s, sf, dx = d->dx;
  distrsum_t *ds;

  /* compute the total number of visits */
  for (tot = 0., i = 0; i < n; i++) tot += d->arr[i].s;

  /* estimate using integral identity */
  for (i = 0; i <= n; i++) {
    /* compute # of visits in a window around i
     * the window is asymmetric at the boundaries */
    j0 = intmax(i - m, 0);
    j1 = intmin(i + m, n);
    for (s0 = s1 = 0., j = j0; j < j1; j++) {
      /* compute the linear phi(x) = (x - x0)/(x1 - x0) - Step(x);
       * here variable `phi' is phi(x) * (x1 - x0)
       * which yields a unit derivative */ 
      phi = (j + .5 - (j < i ? j0 : j1)) * dx;
      ds = d->arr + j;
      s = ds->s;
      sf = ds->sf;
      /* apply the normalization factor, for the radial distribution function */
      if (d->norm) {
        double nm = .5 * (d->norm[j] + d->norm[j + 1]);
        s /= nm;
        sf /= nm;
      }
      s0 += s; /* the histogram part, phi' = 1 */
      s1 += sf * phi;  /* the additive correction */
    }
    den = 1.0 * (j1 - j0) * tot * dx;
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
 * for the volume distribution
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
 * for the volume distribution
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
 *  d->lnrho[i] = \sum_{j < i} d->mf[i] * dx
 * limiting case of distr_ii0 with infinite window size
 * it uses d->mf[i] (which should have been computed already) 
 * it outputs d->rho[i] and d->lnrho[i] */
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

/* fractional integral identity (low level)
 * set rho[], lnrho[], his[]
 * before calling this function
 * first call mfX() to compute the mean force, 
 * and call winX() to compute the window */
INLINE void distr_ii0(distr_t *d, distrwin_t *win)
{
  int i, j, jl, jr, n = d->n;
  double x, delx, den, tot, dx = d->dx;

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
    for (den = 0., delx = 0., j = jl; j <= jr; j++) {
      x = exp(d->lnrho[j] - d->lnrho[i]);
      if (d->norm) x *= d->norm[j]; /* apply the normalization factor */
      den += (j == jl || j == jr) ? x * .5 : x;
      if (d->norm)
        delx += d->norm[j] * ((j == jl || j == jr) ? .5 : 1.);
    }
    den *= tot * dx;
    d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
    if (d->norm == NULL) delx = jr - jl;
    if (tot > 0.) d->his[i] /= tot * dx * delx;
  }
}

/* perform integral identity
 * iitype: 0: Adib-Jazynski, 1: fraction identity, 2: integrate the mean force
 * halfwin: half window size
 *          > 0:  fixed size,
 *            0:  fixed size (determined automatically)
 *           -1:  adaptive window size
 * mfhalfwin: half window size for mean force integration
 *            use zero for simple bin average
 *            positive for unbiased window average 
 *            negative for simple average
 * gam: window size amplification factor
 * mlimit: maximal half window size (-1 to disable)
 * sampmin: minimal number of samples in estimating mean force std. */
INLINE void distr_iiez(distr_t *d, int iitype, int halfwin, int mfhalfwin,
    double gam, int mlimit, double sampmin)
{
  if (iitype == 0) { /* Jarzynski method */
    distr_aj(d, halfwin);
    return;
  } else if (iitype < 0 || iitype > 2) {
    fprintf(stderr, "unsupported method %d\n", iitype);
    return;
  }
  
  /* 1. compute the mean force */
  if (mfhalfwin == 0) {
    distr_mf0(d);
  } else {
    /* positive `mfhalfwin' means an unbiased estimate
     * negative means the */
    distr_mfwin(d, abs(mfhalfwin), mfhalfwin > 0 ? 1 : 0);
  }

  /* 2. compute the distribution */
  if (iitype == 1) {
    if (halfwin < -1) {
      fprintf(stderr, "iiadp no longer supported!\n");
      /* distr_iiadp(d, mlimit, sampmin); */
      return;
    }

    /* computing the window size */
    if (halfwin >= 0) {
      /* uniform windows i.e.
       * if halfwin == 0, estimate from the global mean force fluctuation */
      distr_winfixed(d, d->win, gam, halfwin);
    } else { /* halfwin == -1 */
      /* compute the nonuniform window size */
      distr_mfsig0(d, sampmin);
      distr_winadp0(d, d->win, gam, mlimit);
    }
    /* using the integral identity */
    distr_ii0(d, d->win);
    return;
  }
  
  if (iitype == 0) { /* do simple mean force integration */
    distr_iimf(d);
    return;
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


#define DISTR_KEEPLEFT   0x0040
#define DISTR_KEEPRIGHT  0x0080
#define DISTR_KEEPEDGE   (DISTR_KEEPLEFT | DISTR_KEEPRIGHT)
#define distr_save(d, fn) distr_savex(d, fn, 0)

INLINE int distr_savex(distr_t *d, const char *fn, unsigned flags)
{
  FILE *fp;
  int i, i0 = 0, i1 = d->n, n = d->n;
  distrsum_t *ds = d->arr;

  if (!(flags & DISTR_KEEPEDGE)) {
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

#endif


