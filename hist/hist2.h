#include "util.h"
#include "hist.h"
#ifndef HIST2_H__
#define HIST2_H__

/* two dimensional histograms */
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



INLINE double *hist2getsums_(const double *h, int rows, int n,
    double xmin, double dx, int m, double ymin, double dy, double *sums)
{
  double *xav, *yav, *xxc, *yyc, *xyc, x, y, w;
  int i, j, r;

  xav = sums + rows;
  xxc = sums + rows*2;
  yav = sums + rows*3;
  yyc = sums + rows*4;
  xyc = sums + rows*5;
  for (r = 0; r < rows; r++) {
    sums[r] = xav[r] = xxc[r] = yav[r] = yyc[r] = 0.;
    for (i = 0; i < n; i++) {
      x = xmin + (i+.5)*dx;
      for (j = 0; j < m; j++) {
        y = ymin + (j+.5)*dy;
        w = h[r*n*m + i*m + j];
        sums[r] += w;
        xav[r]  += w*x;
        xxc[r]  += w*x*x;
        yav[r]  += w*y;
        yyc[r]  += w*y*y;
        xyc[r]  += w*x*y;
      }
    }
    if (sums[r] > 1e-5) {
      xav[r] /= sums[r];
      xxc[r] = sqrt(xxc[r]/sums[r] - xav[r]*xav[r]);
      yav[r] /= sums[r];
      yyc[r] = sqrt(yyc[r]/sums[r] - yav[r]*yav[r]);
      xyc[r] = xyc[r]/sums[r] - xav[r]*yav[r];
    }
  }
  return sums;
}


/* write 'rows' 2d n^2 histograms to file */
INLINE int hist2save(const double *h, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  const int version = 1; /* v1 allows different dimension in x and y */
  FILE *fp;
  int i, j, r, imax, imin, jmax, jmin, nm;
  const double *p;
  double *sums, fac, delta;

  if (fn == NULL) fn = "HIST2";
  xfopen(fp, fn, "w", return 1);

  nm = n*m;
  xnew(sums, rows * 6);
  hist2getsums_(h, rows, n, xmin, dx, m, ymin, dy, sums);
  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g %d %g %g | ",
      version, flags, rows, n, xmin, dx, m, ymin, dy);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r]);
  fprintf(fp, " | ");
  for (r = 0; r < rows; r++) /* averages and standard deviations */
    fprintf(fp, "%g %g %g %g %g ", sums[r+rows], sums[r+rows*2],
        sums[r+rows*3], sums[r+rows*4], sums[r+rows*5]);
  fprintf(fp, "| \n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rows; r++) { /* the rth data set */
    p = h + r*nm;

    if (flags & HIST_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge of i */
      for (i = n-1; i >= 0; i--) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imax = i+1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge of i */
      for (i = 0; i < imax; i++) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imin = i;
    }

    if (flags & HIST_KEEPRIGHT2) {
      jmax = m;
    } else { /* trim the right edge of j */
      for (j = m-1; j >= 0; j--) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmax = j+1;
    }

    if (flags & HIST_KEEPLEFT2) {
      jmin = 0;
    } else { /* trim the left edge of j */
      for (j = 0; j < jmax; j++) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmin = j;
    }

    if (fabs(sums[r]) < 1e-6) fac = 1.;
    else fac = 1.0/(sums[r]*dx*dy);

    for (i = imin; i < imax; i++) {
      for (j = jmin; j < jmax; j++) {
        double x, y;
        if ((flags & HIST_NOZEROES) && p[i*m + j] < 1e-16)
          continue;
        x = xmin + (i+delta)*dx;
        y = ymin + (j+delta)*dy;
        fprintf(fp, "%g %g ", x, y);
        if (flags & HIST_KEEPHIST)
          fprintf(fp, "%20.14E ", p[i*m+j]);
        fprintf(fp, "%20.14E %d\n", p[i*m+j]*fac, r);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n#\n");
  }
  fclose(fp);
  if (flags & HIST_VERBOSE) {
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f xav: %10.4f(%10.4f) yav: %10.4f(%10.4f)\n",
          r, sums[r], sums[r+rows], sums[r+rows*2], sums[r+rows*3], sums[r+rows*4]);
  }
  free(sums);
  return 0;
}

INLINE int hist2load(double *hist, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  FILE *fp;
  static char s[40960] = "", *p;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, j, r, r1, nm, nlin = 0;
  unsigned fflags;
  double x, y, g, g2, fac, delta, *arr, *sums = NULL;
  double xmin1, dx1, ymin1, dy1;

  xfopen(fp, fn, "r", return -1);

  nm = n*m;
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;
  if (9 != sscanf(s, " # %d 0x %X | %d%d%lf%lf%d%lf%lf | %n", &ver, &fflags, &r,
        &i, &xmin1, &dx1, &j, &ymin1, &dy1, &next)
      || i < n || j < m || r != rows
      || fabs(dx1 - dx) > 1e-5 || fabs(dy1 - dy) > 1e-5 ) {
    fprintf(stderr, "Error: bins %d, %d; %d, %d; ng %d, %d; dx %g, %g; dy %g, %g\n",
        i, n, j, m, r, rows, dx1, dx, dy1, dy);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST_KEEPHIST);
  /* scan sums */
  xnew(sums, rows);
  for (p = s+next, r = 0; r < rows; r++) {
    if (1 != sscanf(p, "%lf%n", sums + r, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;
  for (r = 0; r < rows; r++) {
    if (4 != sscanf(p, "%lf%lf%lf%lf%n", &x, &y, &g, &g2, &next)) {
      fprintf(stderr, "cannot read average/stddev from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;

  if (!add) { /* clear histogram */
    for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  }

  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    arr = hist + r*nm;
    fac = sums[r]*(dx*dx);
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; isspace((unsigned char)(*p)) && p >= s; p--)
        *p = '\0'; /* trim ending */
      if (s[0] == '#') break;
      if (s[0] == '\0') continue;

      if (hashist) {
        if (5 != sscanf(s, "%lf%lf%lf%lf%d", &x, &y, &g, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if (r1 != r) {
        fprintf(stderr, "wrong column index %d vs. %d on line %d\n",
          r1, r, nlin);
        goto EXIT;
      }
      i = (int)((x - xmin)/dx - delta + .5);
      if (i < 0 || i >= n) {
        fprintf(stderr, "cannot find index for x = %g\n", x);
        goto EXIT;
      }
      j = (int)((y - ymin)/dy - delta + .5);
      if (j < 0 || j >= m) {
        fprintf(stderr, "cannot find index for y = %g\n", y);
        return -1;
      }
      if (!hashist) {
        g = g2*fac;
      }
      if (add) arr[i*m+j] += g;
      else arr[i*m+j] = g;
    }
  }
  if (verbose) fprintf(stderr, "%s loaded successfully\n", fn);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  if (sums) free(sums);
  for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  return -1;
}

/* add (xarr[skip*r], yarr[skip*r]) of weight w, into histogram h
 * return number of success */
INLINE int hist2add(const double *xarr, const double *yarr, int skip,
    double w, double *h, int rows,
    int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags)
{
  int r, ix, iy, good = 0, verbose = flags & HIST_VERBOSE;
  double x, y;

  for (r = 0; r < rows; r++) {
    x = xarr[skip*r];
    y = yarr[skip*r];
    if (x < xmin || y < ymin) {
      if (verbose)
        fprintf(stderr, "histadd underflows %d: %g < %g or %g < %g\n",
          r, x, xmin, y, ymin);
      continue;
    }
    ix = (int)((x - xmin)/dx);
    iy = (int)((y - ymin)/dy);
    if (ix >= n || iy >= m) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g > %g or %g > %g\n",
            r, x, xmin + dx*n, y, ymin + dy*m);
      continue;
    }
    h[r*n*m + ix*m + iy] += w;
    good++;
  }
  return good;
}


#endif /* ZCOM_HIST2__ */

