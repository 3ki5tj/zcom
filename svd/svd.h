#include "util.h"
#ifndef SVD_H__
#define SVD_H__



/* singular value decomposition of mxn matrix `a'
 * a[m*n] (or u[m*n] on return), w[n], v[n*n] */
INLINE int svd(real *a, real *w, real *v, int m, int n)
{
  int flag, i, it, j, jj, k, l, nm;
  real c, f, h, s, x, y, z;
  real anorm = 0, g, scl;
  real *rv1;

  die_if (m < n, "ERROR: m %d < n %d\n", m, n);
  xnew(rv1, n);

  /* Householder reduction to bidiagonal form */
  for (g = s = scl = 0., i = 0; i < n; i++) {
    /* left-hand reduction */
    l = i + 1;
    rv1[i] = scl * g;
    g = s = scl = 0.0;
    if (i < m) {
      for (k = i; k < m; k++)
        scl += (real) fabs(a[k*n+i]);
      if (scl > 0.) {
        for (k = i; k < m; k++) {
          a[k*n+i] = x = a[k*n+i]/scl;
          s += x*x;
        }
        f = a[i*n+i];
        g = (real) ((f > 0.) ? -sqrt(s) : sqrt(s));
        h = f * g - s;
        a[i*n+i] = f - g;
        if (i != n - 1) {
          for (j = l; j < n; j++) {
            for (s = 0.0, k = i; k < m; k++)
              s += a[k*n+i] * a[k*n+j];
            f = s / h;
            for (k = i; k < m; k++)
              a[k*n+j] += f * a[k*n+i];
          }
        }
        for (k = i; k < m; k++)
          a[k*n+i] = a[k*n+i]*scl;
      }
    }
    w[i] = scl*g;

    /* right-hand reduction */
    g = s = scl = 0.0;
    if (i < m && i != n - 1) {
      for (k = l; k < n; k++)
        scl += (real) fabs(a[i*n+k]);
      if (scl > 0.) {
        for (k = l; k < n; k++) {
          a[i*n+k] = x = a[i*n+k]/scl;
          s += x*x;
        }
        f = a[i*n+l];
        g = (real) ((f > 0) ? -sqrt(s) : sqrt(s));
        h = f * g - s;
        a[i*n+l] = f - g;
        for (k = l; k < n; k++)
          rv1[k] = a[i*n+k] / h;
        if (i != m - 1) {
          for (j = l; j < m; j++) {
            for (s = 0.0, k = l; k < n; k++)
              s += a[j*n+k] * a[i*n+k];
            for (k = l; k < n; k++)
              a[j*n+k] += s * rv1[k];
          }
        }
        for (k = l; k < n; k++)
          a[i*n+k] *= scl;
      }
    }
    x = (real) fabs(w[i]);
    x += (real) fabs(rv1[i]);
    if (x > anorm) anorm = x;
  }

  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--) {
    if (i < n - 1) {
      if (g != 0) {
        for (j = l; j < n; j++)
          v[j*n+i] = ((a[i*n+j] / a[i*n+l]) / g);
        /* real division to avoid underflow */
        for (j = l; j < n; j++) {
          for (s = 0.0, k = l; k < n; k++)
            s += (a[i*n+k] * v[k*n+j]);
          for (k = l; k < n; k++)
            v[k*n+j] += (s * v[k*n+i]);
        }
      }
      for (j = l; j < n; j++)
        v[i*n+j] = v[j*n+i] = 0.0;
    }
    v[i*n+i] = 1.0;
    g = rv1[i];
    l = i;
  }

  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    if (i < n - 1)
      for (j = l; j < n; j++) a[i*n+j] = 0;
    if (g != 0) {
      g = 1 / g;
      if (i != n - 1) {
        for (j = l; j < n; j++) {
          for (s = 0.0, k = l; k < m; k++)
            s += (a[k*n+i] * a[k*n+j]);
          f = s/a[i*n+i]*g;
          for (k = i; k < m; k++)
            a[k*n+j] += f*a[k*n+i];
        }
      }
      for (j = i; j < m; j++)
        a[j*n+i] = a[j*n+i]*g;
    } else {
      for (j = i; j < m; j++) a[j*n+i] = 0.0;
    }
    a[i*n+i] += 1.;
  }

  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--) { /* loop over singular values */
    for (it = 0; it < 200; it++) { /* loop over allowed iterations */
      flag = 1;
      for (l = k; l >= 0; l--) { /* test for splitting */
        nm = l - 1;
        if (fabs(rv1[l]) + anorm == anorm) {
          flag = 0;
          break;
        }
        if (fabs(w[nm]) + anorm == anorm)
          break;
      }
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          if (fabs(f) + anorm == anorm) continue;
          g = w[i];
          h = (real) dblhypot(f, g);
          w[i] = h;
          h = 1 / h;
          c = g * h;
          s = (- f * h);
          for (j = 0; j < m; j++) {
            y = a[j*n+nm];
            z = a[j*n+i];
            a[j*n+nm] = y * c + z * s;
            a[j*n+i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) { /* convergence */
        if (z < 0.0) { /* flip sign of w */
          w[k] = -z;
          for (j = 0; j < n; j++)
            v[j*n+k] = -v[j*n+k];
        }
        break;
      }
      if (it >= 200) {
        free(rv1);
        fprintf(stderr, "svd: failed to converge\n");
        return -1;
      }

      /* shift from bottom 2 x 2 minor */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
      g = (real) dblhypot(f, 1.0);
      if (f < 0.) g = -g;
      f = ((x - z) * (x + z) + h * (y/(f + g) - h)) / x;

      /* next QR transformation */
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = (real) dblhypot(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++) {
          x = v[jj*n+j];
          z = v[jj*n+i];
          v[jj*n+j] = x * c + z * s;
          v[jj*n+i] = z * c - x * s;
        }
        w[j] = z = (real) dblhypot(f, h);
        if (z > 0.) { c = f/z; s = h/z; }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < m; jj++) {
          y = a[jj*n+j];
          z = a[jj*n+i];
          a[jj*n+j] = y * c + z * s;
          a[jj*n+i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free(rv1);
  return 0;
}



INLINE int svdback(real *u, real *w, real *v, int m, int n, real *x, real *b)
{
  int i, j;
  real *b1, y;

  xnew(b1, n);
  for (i = 0; i < n; i++) {
    if (w[i] <= 0.) { b1[i] = 0.; continue; }
    for (y = 0, j = 0; j < m; j++)
      y += u[j*n+i]*b[j];
    b1[i] = y/w[i];
  }
  for (i = 0; i < n; i++) {
    for (y = 0., j = 0; j < n; j++)
      y += v[i*n+j]*b1[j];
    x[i] = y;
  }
  free(b1);
  return 0;
}



INLINE int svdsolve(real *a, real *x, real *b, int n, real rerr)
{
  int i;
  real *u, *v, *w, wmax, wmin;

  xnew(w, n); xnew(u, n*n); xnew(v, n*n);
  for (i = 0; i < n*n; i++) u[i] = a[i];
  svd(u, w, v, n, n);
  for (wmax = 0., i = 0; i < n; i++)
    if (w[i] > wmax) wmax = w[i];
  for (wmin = wmax*rerr, i = 0; i < n; i++)
    if (w[i] < wmin) w[i] = wmin;
  for (i = 0; i < n; i++) printf("%g  ", w[i]);
  printf("\n");
  svdback(u, w, v, n, n, x, b);
  free(u); free(v); free(w);
  return 0;
}

#endif

