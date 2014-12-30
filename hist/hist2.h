#include "util.h"
#include "hist.h"
#ifndef HIST2_H__
#define HIST2_H__



/* two dimensional histograms */



#define HIST2_VERBOSE    0x0001
#define HIST2_ADDAHALF   0x0010
#define HIST2_NOZEROES   0x0020
#define HIST2_KEEPLEFT   0x0040
#define HIST2_KEEPRIGHT  0x0080
#define HIST2_KEEPLEFT2  0x0040
#define HIST2_KEEPRIGHT2 0x0080
#define HIST2_KEEPEDGE   (HIST2_KEEPLEFT | HIST2_KEEPRIGHT | HIST2_KEEPLEFT2 | HIST2_KEEPRIGHT2)
#define HIST2_KEEPHIST   0x0100
#define HIST2_OVERALL    0x0200
#define HIST2_ADDITION   0x1000



INLINE void hist2getsums_(const double *h, int n,
    double xmin, double dx, int m, double ymin, double dy, double *s)
{
  double x, y, w;
  int i, j;

  s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0;
  for (i = 0; i < n; i++) {
    x = xmin + (i+.5)*dx;
    for (j = 0; j < m; j++) {
      y = ymin + (j+.5)*dy;
      w = h[i*m + j];
      s[0]  += w;
      s[1]  += w * x;
      s[2]  += w * y;
      s[3]  += w * x * x;
      s[4]  += w * x * y;
      s[5]  += w * y * y;
    }
  }
  if ( s[0] > 0 ) {
    s[1] /= s[0];
    s[2] /= s[0];
    s[3]  = s[3] / s[0] - s[1] * s[1];
    s[4]  = s[4] / s[0] - s[1] * s[2];
    s[5]  = s[5] / s[0] - s[2] * s[2];
  }
}



/* write rows n x m histograms to file */
INLINE int hist2save(const double *h, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  const int version = 1; /* v1 allows different dimension in x and y */
  FILE *fp;
  int i, j, r, imax, imin, jmax, jmin, nm;
  const double *p;
  double (*sums)[6], fac, delta;

  if (fn == NULL) fn = "HIST2";
  xfopen(fp, fn, "w", return 1);

  nm = n * m;
  xnew(sums, rows);
  for ( r = 0; r < rows; r++ )
    hist2getsums_(h + nm, n, xmin, dx, m, ymin, dy, sums[r]);
  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g %d %g %g | ",
      version, flags, rows, n, xmin, dx, m, ymin, dy);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r][0]);
  fprintf(fp, " | ");
  for (r = 0; r < rows; r++) /* averages and covariance */
    fprintf(fp, "%g %g %g %g %g ", sums[r][1], sums[r][2],
        sums[r][3], sums[r][4], sums[r][5]);
  fprintf(fp, "| \n");

  delta = (flags & HIST2_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rows; r++) { /* the rth data set */
    p = h + r*nm;

    if (flags & HIST2_KEEPRIGHT) {
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

    if (flags & HIST2_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge of i */
      for (i = 0; i < imax; i++) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imin = i;
    }

    if (flags & HIST2_KEEPRIGHT2) {
      jmax = m;
    } else { /* trim the right edge of j */
      for (j = m-1; j >= 0; j--) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmax = j+1;
    }

    if (flags & HIST2_KEEPLEFT2) {
      jmin = 0;
    } else { /* trim the left edge of j */
      for (j = 0; j < jmax; j++) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmin = j;
    }

    if ( fabs(sums[r][0]) < 1e-6 ) fac = 1.;
    else fac = 1.0 / (sums[r][0] * dx * dy);

    for (i = imin; i < imax; i++) {
      for (j = jmin; j < jmax; j++) {
        double x, y;
        if ((flags & HIST2_NOZEROES) && p[i*m + j] < 1e-16)
          continue;
        x = xmin + (i+delta)*dx;
        y = ymin + (j+delta)*dy;
        fprintf(fp, "%g %g ", x, y);
        if (flags & HIST2_KEEPHIST)
          fprintf(fp, "%20.14E ", p[i*m+j]);
        fprintf(fp, "%20.14E %d\n", p[i*m+j]*fac, r);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n#\n");
  }
  fclose(fp);
  if (flags & HIST2_VERBOSE) {
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f xav: %10.4f(%10.4f) yav: %10.4f(%10.4f)\n",
          r, sums[r][0], sums[r][1], sums[r][2], sums[r][3], sums[r][4]);
  }
  free(sums);
  return 0;
}



INLINE int hist2load(double *hist, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  FILE *fp;
  int verbose = (flags & HIST2_VERBOSE);
  int add = (flags & HIST2_ADDITION);
  int ver, next, hashist;
  int i, j, r, r1, nm, nlin = 0;
  unsigned fflags;
  double x, y, xx, xy, yy, g, g2, fac, delta, *arr, *sums = NULL;
  double xmin1, dx1, ymin1, dy1;
  char s[40960] = "", *p;

  xfopen(fp, fn, "r", return -1);

  nm = n * m;
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
  delta   = ((fflags & HIST2_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST2_KEEPHIST);
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
    if (5 != sscanf(p, "%lf%lf%lf%lf%lf%n", &x, &y, &xx, &xy, &yy, &next)) {
      fprintf(stderr, "cannot read ave./cov. from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;

  if ( !add ) { /* clear histogram */
    for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  }

  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    arr = hist + r*nm;
    fac = sums[r] * dx * dy;
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; p >= s && isspace((unsigned char)(*p)); p--)
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
  fclose(fp);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  if (sums) free(sums);
  for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  fclose(fp);
  return -1;
}



/* add (x, y) into the histogram h with weight w
 * return the number of successful rows */
INLINE int hist2add1(double x, double y, double w, double *h,
    int n, double xmin, double dx,
    int m, double ymin, double dy, int verbose)
{
  int ix, iy;

  if ( x < xmin ) {
    if ( verbose )
      fprintf(stderr, "hist2add underflows: x %g < %g\n",
          x, xmin);
    return -1;
  }
  if ( y < xmin ) {
    if ( verbose )
      fprintf(stderr, "hist2add underflows: y %g < %g\n",
          y, ymin);
    return -1;
  }
  if ( (ix = (int)((x - xmin)/dx)) >= n ) {
    if (verbose)
      fprintf(stderr, "hist2add overflows: x %g > %g\n",
          x, xmin + dx*n);
    return -1;
  }
  if ( (iy = (int)((y - ymin)/dy)) >= m ) {
    if (verbose)
      fprintf(stderr, "hist2add overflows: y %g > %g\n",
          y, ymin + dy*m);
    return -1;
  }
  h[ix*m + iy] += w;
  return 0;
}



/* add (xarr[stride * r], yarr[stride * r]) of weight w
 * into the rth row of histogram h, r = 0..rows - 1
 * return the number of successful rows */
INLINE int hist2add(const double *xarr, const double *yarr, int stride,
    double w, double *h, int rows,
    int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags)
{
  int r, ix, iy, good = 0, verbose = flags & HIST2_VERBOSE;
  double x, y;

  for (r = 0; r < rows; r++) {
    if ( (x = xarr[stride * r]) < xmin ) {
      if ( verbose )
        fprintf(stderr, "hist2add underflows, row %d: x %g < %g\n",
            r, x, xmin);
      continue;
    }
    if ( (y = yarr[stride * r]) < xmin ) {
      if ( verbose )
        fprintf(stderr, "hist2add underflows, row %d: y %g < %g\n",
            r, y, ymin);
      continue;
    }
    if ( (ix = (int)((x - xmin)/dx)) >= n ) {
      if (verbose)
        fprintf(stderr, "hist2add overflows, row %d: x %g > %g\n",
            r, x, xmin + dx*n);
      continue;
    }
    if ( (iy = (int)((y - ymin)/dy)) >= m ) {
      if (verbose)
        fprintf(stderr, "hist2add overflows, row %d: y %g > %g\n",
            r, y, ymin + dy*m);
      continue;
    }
    h[r*n*m + ix*m + iy] += w;
    good++;
  }
  return good;
}



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
  free(hs2->arr);
  free(hs2);
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



#define hs2_add1ez(hs, x, y, flags) hs2_add1(hs, 0, x, y, 1.0, flags)

INLINE int hs2_add1(hist2_t *hs, int r, double x, double y,
    double w, unsigned flags)
{
  hs2_check(hs);
  return hist2add1(x, y, w, hs->arr+r*hs->n*hs->m,
      hs->n, hs->xmin, hs->dx, hs->m, hs->ymin, hs->dy, flags & HIST2_VERBOSE);
}



INLINE int hs2_add(hist2_t *hs, const double *x, const double *y, int stride,
    double w, unsigned flags)
{
  int r, good = 0;
  hs2_check(hs);
  for ( r = 0; r < hs->rows; r++ )
    good += (hs2_add1(hs, r, x[stride*r], y[stride*r], w, flags) == 0);
  return good;
}



#endif /* ZCOM_HIST2__ */

