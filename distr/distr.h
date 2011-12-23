#include "def.h"
#include "util.h"
#include "rc.h"
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
  printf("mean force standard deviation is %g\n", sig);

  m = intmax(1, (int)( 0.5/(sig * d->dx) + .5 ));
  distr_winfix(d, m);
}

/* estimate local standard deviation of the mean force
 * average over local window of 2 m */
INLINE void distr_mfsig1(distr_t *d, double sampmin)
{
  int i, m, n = d->n;
  double x, sig, his, *asig, *hsum;

  xnew(asig, n + 1);
  xnew(hsum, n + 1);
  for (asig[0] = hsum[0] = 0, i = 0; i < n; i++) {
    distrsum_t *ds = d->arr + i;
    hsum[i+1] = hsum[i] + ds->s;
    x = ds->sf2 * ds->s - ds->sf * ds->sf;
    asig[i+1] = asig[i] + ((ds->s > 0. && x > 0.) ? sqrt(x) : 0.);
  }
  for (i = 0; i < n; i++) {
    /* start from the minimal window size, expand until included sampmin data points */
    for (m = 1; m < n; m++) {
      int j0 = intmax(0, i-m), j1 = intmin(i + m + 1, n);
      his = hsum[j1] - hsum[j0];
      sig = asig[j1] - asig[j0];
      d->mfsig[i] = his > 0. ? sig/his : 0.;
      if (his > sampmin || (j0 == 0 && j1 == n)) break;
    }
  }
  d->mfsig[n] = d->mfsig[n-1];
  free(asig);
  free(hsum);
} 

/* make adaptive window
 * requires d->rho and d->mfsig */
INLINE void distr_winadp(distr_t *d)
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

  for (i = 0; i <= n; i++) {
    int k = 0, jl = i, jr = i+1, el = 0, er = 0;

    for (k = 0; el == 0 || er == 0; k++) {
      if (k % 2 == 0) { /* do the left */
        if (el || jl <= 1) {
          el = 1;
          continue;
        }
        if (hhis[jl - 1] / (hsum[jr] - hsum[jl] + 1e-10) < d->mfsig[jl - 1] * dx) {
          el = 1;
        } else {
          jl--;
        }
      } else { /* do the right */
        if (er || jr >= n) {
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
  int i, n = d->n;
  double max, tot, lntot, dx = d->dx;

  /* construct ln(rho) from the mean force */
  for (max = d->lnrho[0] = 0, i = 0; i < n; i++) {
    d->lnrho[i + 1] = d->lnrho[i] + d->mf[i] * dx;
    if (d->lnrho[i + 1] > max) max = d->lnrho[i + 1];
  }

  for (tot = 0, i = 0; i <= n; i++)
    tot += d->rho[i] = exp(d->lnrho[i] -= max);
  tot *= dx;
  lntot = log(tot);

  for (i = 0; i <= n; i++) {
    d->lnrho[i] -= lntot;
    d->rho[i] /= tot;
  }
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
      if (d->norm) x *= d->norm[j];
      den += (j == jl || j == jr) ? x * .5 : x;
    }
    den *= tot * dx;
    d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
    d->his[i] /= tot * dx * (jr - jl);
  }
}


INLINE void distr_iiadp(distr_t *d, double sigsampmin)
{
  distr_iimf(d); /* tentative distribution from integrating the mean force */
  distr_mfsig1(d, sigsampmin);
  distr_ii0(d, -1); /* use variable window */
}

/* perform integral identity
 * iitype: 0: Adib-Jazynski, 1: fraction identity, 2: integrate the mean force
 * halfwin: half window size. > 0: fixed size, 0: fixed size (determined
 *          automatically), -1: adaptive window size
 * mfhalfwin: half window size for mean force integration
 * sampmin: minimal number of samples in estimating mean force std. */
INLINE void distr_iiez(distr_t *d, int iitype, int halfwin, int mfhalfwin,
    double sampmin)
{
  if (iitype == 0) {
    distr_aj(d, halfwin);
  } else if (iitype == 1 || iitype == 2) {
    if (mfhalfwin == 0) {
      distr_mf0(d);
    } else if (mfhalfwin > 0) {
      distr_mfii(d, mfhalfwin);
    } else if (mfhalfwin < 0) {
      distr_mfav(d, -mfhalfwin);
    }

    if (iitype == 1) {
      if (halfwin >= 0)
        distr_ii0(d, halfwin);
      else
        distr_iiadp(d, sampmin);
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
  xnew(d->jl, d->n + 1);
  xnew(d->jr, d->n + 1);
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
  free(d->jl);
  free(d->jr);
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
      d->xmin + d->dx * d->jl[i], d->xmin + d->dx * d->jr[i], 
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
      x = rc_smul(fts->rf[kl], sk * cl * dx);
      x = rc_add(x, rc_smul(fts->rg[kl], ck * sl * dy));
      x = rc_mul(x, rc_make(sk*cl - ck*sl, ck*cl + sk*sl));
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

typedef struct {
  int n, m;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
  distr2dsum_t *arr;
  double *rho, *lnrho, *his;
  double *mf, *mg;
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

  n1m1 = (d->n + 1) * (d->m + 1);
  xnew(d->arr, n1m1);
  xnew(d->rho, n1m1);
  xnew(d->lnrho, n1m1);
  xnew(d->his, n1m1);
  xnew(d->mf, n1m1);
  xnew(d->mg, n1m1);
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

/* compute the mean force */
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

/* compute lnrho from least-square fitting of the mean force */
INLINE void distr2d_calclnrho(distr2d_t *d)
{
  ftsolver2d_t *fts = ftsolver2d_open(d->n, (real) d->dx, d->m, (real) d->dy);
  ftsolver2d_solve(fts, d->lnrho, d->mf, d->mg);
  ftsolver2d_close(fts);
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
