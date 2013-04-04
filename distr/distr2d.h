#ifndef DISTR2D_H__
#define DISTR2D_H__

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
  int i, j, n = d->n, m = d->m, m1 = d->m + 1, ww;

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
        if (ms > smin) break; 
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
      fprintf(fp, "%g %g %.14e %.14e %.14e %.14e %.14e %.14e %g %g %g %d\n",
        d->xmin + i * d->dx, d->ymin + j * d->dy, 
        sm, f, g, f2, g2, fg,
        d->rho[id], d->lnrho[id], d->his[id], d->win[id].ir - d->win[id].il);
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

