#include "lu.h"
#ifndef SAVGOL_C__
#define SAVGOL_C__

#include "savgol.h"

/* compute 1d Savitzky-Golay coefficients 
 * der == 0 for function itself, 1 for first-order derivative */
INLINE double *savgol(int w, int ord, int der, int h, int verbose)
{
  int i, i0, i1, ox, oy, nop, orm, npt;
  double x, xk, y;
  double *mmt, *b, *mat, *x2m, *c;

  nop = ord+1;
  orm = 2*ord+1;
  npt = h ? (2*w) : (2*w+1);
  i0 = -w;
  i1 = h ? w : (w+1);
  if ((c = calloc(npt, sizeof(double))) == NULL) 
    return NULL;
  if ((b = calloc(nop, sizeof(double))) == NULL) {
    free(c);
    return NULL;
  }
  if ((mat = calloc(nop*nop, sizeof(double))) == NULL) {
    free(b); free(c);
    return NULL;
  }
  if ((mmt = calloc(orm, sizeof(double))) == NULL) {
    free(mat); free(b); free(c);
    return NULL;
  }
  if ((x2m = calloc(nop*npt, sizeof(double))) == NULL) {
    free(mmt); free(mat); free(b); free(c);
    return NULL;
  }
  if (der > ord) {
    fprintf(stderr, "no %dth order of derivative, order = %d\n", der, ord);
    return NULL;
  }

  for (i = 0; i < orm; i++) mmt[i] = 0.;
  for (i = i0; i < i1; i++) {
    x = h ? (i + .5) : i;
    /* mmt[k] = < x^k > */
    for (xk = 1., ox = 0; ox < orm; ox++, xk *= x)
      mmt[ox] += xk;
    /* x2m[k*npt + x] = x^k */
    for (xk = 1., ox = 0; ox <= ord; ox++, xk *= x)
      x2m[ox*npt + (i - i0)] = xk;
  }

  /* install matrix from moments */
  for (ox = 0; ox < nop; ox++)
    for (oy = 0; oy < nop; oy++)
      mat[ox*nop + oy] = mmt[ox+oy];
 
  /* we approximate y(x) = a0 + a1 x + a2 x^2 + ...
   * mat.a = b = x2m.y, or a = mat^(-1).b
   * since mat is symmetrical, rows == columns,
   * we thus extract the first row by solving b = {1, 0, 0, ...} */ 
  for (i = 0; i < nop; i++) b[i] = 0;
  b[der] = 1.;
  i = lusolve(mat, b, nop);
  if (i != 0) {
    fprintf(stderr, "unable to inverse matrix!\n");
    return NULL;
  }
  /* c = mat^(-1).x2m */
  for (i = 0; i < npt; i++) {
    for (y = 0, ox = 0; ox < nop; ox++)
      y += b[ox]*x2m[ox*npt + i];
    c[i] = y;
  }
  free(x2m);
  free(mmt);
  free(mat);
  free(b);
  if (verbose) {
    for (i = 0; i < npt; i++)
      printf("%g\t", c[i]);
    printf("\n");
  }
  return c;
}

/* compute 2d Savitzky-Golay coefficients
 * h means if it is a histogram */
INLINE double *savgol2d(int iw, int jw, int ord, int h, int verbose)
{
  int i, j, i0, i1, j0, j1, id, nop, orm, npt;
  int io, iq, ox, oy, o1, o2, o3, o4;
  double w, x, y, xk, xyk;
  double *mmt, *b, *mat, *x2m, *c;

  nop = (ord+1)*(ord+2)/2;
  orm = 2*ord+1;
  i0 = -iw;
  j0 = -jw;
  if (h) { /* for histogram */
    npt = (2*iw)*(2*jw);
    i1 = iw;
    j1 = jw;
  } else {
    npt = (2*iw+1)*(2*jw+1);
    i1 = iw + 1;
    j1 = jw + 1;
  }
  if ((c = calloc(npt, sizeof(double))) == NULL) 
    return NULL;
  if ((b = calloc(nop, sizeof(double))) == NULL) {
    free(c);
    return NULL;
  }
  if ((mat = calloc(nop*nop, sizeof(double))) == NULL) {
    free(b); free(c);
    return NULL;
  }
  if ((mmt = calloc(orm*orm, sizeof(double))) == NULL) {
    free(mat); free(b); free(c);
    return NULL;
  }
  if ((x2m = calloc(nop*npt, sizeof(double))) == NULL) {
    free(mmt); free(mat); free(b); free(c);
    return NULL;
  }

  for (i = 0; i < orm*orm; i++) mmt[i] = 0.;
  for (i = i0; i < i1; i++) {
    x = h ? (i + .5) : i;
    for (j = j0; j < j1; j++) {
      y = h ? (j + .5) : j;
      w = 1.;
      /* moment matrix */
      xk = w;
      for (ox = 0; ox < orm; ox++) {
        xyk = xk;
        for (oy = 0; oy < orm-ox; oy++) {
          mmt[ox*orm + oy] += xyk;
          xyk *= y;
        }
        xk *= x;
      }
      /* position to z-moment matrix */
      id = (i - i0)*(j1 - j0) + (j - j0);
      for (io = 0, o1 = 0; o1 <= ord; o1++) {
        for (o2 = 0; o2 <= o1; o2++, io++) {
          xyk = w;
          for (ox = 0; ox < o1 - o2; ox++) xyk *= x;
          for (oy = 0; oy < o2; oy++) xyk *= y; 
          x2m[io*npt + id] = xyk;
        }
      }
    }
  }

  /* install matrix from moments */ 
  for (io = 0, o1 = 0; o1 <= ord; o1++)
  for (o2 = 0; o2 <= o1; o2++, io++) {
    /* x^(o1-o2) y^o2 */
    for (iq = 0, o3 = 0; o3 <= ord; o3++)
    for (o4 = 0; o4 <= o3; o4++, iq++) {
      /* x^(o3-o4) y^o4 */
      ox = o3 - o4 + o1 - o2;
      oy = o4 + o2;
      mat[io*nop + iq] = mmt[ox*orm + oy];
    }
  } 
  
  for (i = 0; i < nop; i++) b[i] = 0.;
  b[0] = 1.;
  i = lusolve(mat, b, nop);
  if (i != 0) {
    fprintf(stderr, "unable to inverse matrix!\n");
    return NULL;
  }
  for (i = 0; i < npt; i++) {
    for (y = 0, io = 0; io < nop; io++)
      y += b[io]*x2m[io*npt + i];
    c[i] = y;
  }
  free(x2m);
  free(mmt);
  free(mat);
  free(b);
  if (verbose) {
    for (i = i0; i < i1; i++) {
      for (j = j0; j < j1; j++) {
        printf("%7.4f ", c[(i-i0)*(j1 - j0)+(j-j0)]);
      }
      printf("\n");
    }
  }
  return c;
}

#endif

