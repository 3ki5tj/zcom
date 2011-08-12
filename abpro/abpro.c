#include "util.c"

#ifndef ABPRO_C__
#define ABPRO_C__
/* AB beads protein models */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "abpro.h"

/* initialization
 * seqid: 8: 34, 9: 55, 10: 89*/
abpro_t *ab_open(int seqid, int d, int model)
{
  abpro_t *ab;
  int i, nd;
  double x;

  die_if (d == 2 && model != 1, "2d only for model 1");
  die_if (seqid < 0, "bad seqid %d\n", seqid);

  xnew(ab, 1);
  ab->d = d;

  ab->model = model;
  if (model == 1) {
    ab->clj[0][0] = 1; ab->clj[1][1] = .5f;
    ab->clj[0][1] = ab->clj[1][0] = -.5f;
    ab->sla = ab->slb = 24;
  } else {
    ab->clj[0][0] = 1;
    ab->clj[0][1] = ab->clj[1][0] = ab->clj[1][1] = .5f;
    ab->sla = ab->slb = 24;
  }

  /* determine # of atoms */
  x = pow(.5*(sqrt(5.) + 1), seqid + 1);
  i = (seqid % 2) ? (-1) : 1;
  ab->n = (int)( (x + i/x)/sqrt(5.) + .5);

  /* construct sequence */
  xnew(ab->type, ab->n);
  if (seqid < 2) {
    ab->type[0] = seqid;
  } else {
    int *s[2], sl[2], who, j;
    xnew(s[0], ab->n);
    xnew(s[1], ab->n);
    s[0][0] = 0; sl[0] = 1;
    s[1][0] = 1; sl[1] = 1;
    for (who = 0, i = 2; i <= seqid; i++, who = !who) {
      /* s[who] += s[!who]; */
      die_if (sl[0] + sl[1] > ab->n, "combined length > %d\n", ab->n);
      for (j = 0; j < sl[!who]; j++)
        s[who][sl[who] + j] = s[!who][j];
      sl[who] += sl[!who];
    }
    for (who = !who, j = 0; j < ab->n; j++) {
     ab->type[j] = s[who][j];
    }
    free(s[0]);
    free(s[1]);
  }
  printf("%d: ", ab->n);
  for (i = 0; i < ab->n; i++) printf("%c", ab->type[i]+'A'); printf("\n");

  nd = ab->n * ab->d;
  xnew(ab->x, nd);
  xnew(ab->x1, nd);
  xnew(ab->dx, nd);
  xnew(ab->dx1, nd);
  xnew(ab->v, nd);
  xnew(ab->f, nd);
  xnew(ab->lmx, nd);
  xnew(ab->xmin, nd);
  return ab;
}

/* close ab */
void ab_close(abpro_t *ab)
{
  if (ab) {
    ab_milcshake(ab, NULL, NULL, NULL, 0, 0, 0);
    free(ab->type);
    free(ab->x);
    free(ab->x1);
    free(ab->dx);
    free(ab->dx1);
    free(ab->v);
    free(ab->f);
    free(ab->lmx);
    free(ab->xmin);
    free(ab);
  }
}

/* check connectivity */
int ab_checkconn(abpro_t *ab, const real *x)
{
  int i, d = ab->d;
  double r;

  for (i = 0; i < ab->n-1; i++) {
    if (d == 3) {
      r = rv3_dist(x + i*3, x + (i+1)*3);
    } else {
      r = rv2_dist(x + i*2, x + (i+1)*2);
    }
    if (fabs(r-1) > 1e-3) {
      fprintf(stderr, "link (%d,%d) is broken.\n", i, i+1);
      return 1;
    }
  }
  return 0;
}

/* shift the center of mass to zero */
void ab_shiftcom(abpro_t *ab, real *x)
{
  int i, j, d = ab->d, n = ab->n;
  real rc;

  for (j = 0; j < d; j++) {
    rc = 0;
    for (i = 0; i < n; i++) rc += x[i*d+j];
    rc /= n;
    for (i = 0; i < n; i++) x[i*d+j] -= rc;
  }
}

/* write position file (which may include velocity) */
int ab_writepos(abpro_t *ab, const real *x, const real *v, const char *fname)
{
  FILE *fp;
  int i, j, d = ab->d, n = ab->n;

  if (fname == NULL) fname = "ab.pos";
  if ((fp = fopen(fname, "w")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fname);
    return 1;
  }

  fprintf(fp, "# %d %d %d %d\n", d, n, ab->model, (v != NULL));
  for (i = 0; i < n; i++) {
    for (j = 0; j < d; j++) fprintf(fp, "%16.14f ", x[i*d+j]);
    if (v)
      for (j = 0; j < d; j++) fprintf(fp, "%16.14f ", v[i*d+j]);
    fprintf(fp, "%d ", ab->type[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}


/* read position file (which may include velocity) */
int ab_readpos(abpro_t *ab, real *x, real *v, const char *fname)
{
  char s[1024], *p;
  FILE *fp;
  int i, j, hasv = 0, next, d = ab->d, n = ab->n;
  const char *fmt;
  real vtmp[3], *vi;

  if (fname == NULL) fname = "ab.pos";
  if ((fp = fopen(fname, "r")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fname);
    return 1;
  }

  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "Warning: %s has no information line\n", fname);
    rewind(fp);
  } else {
    if (3 != sscanf(s+1, "%d%d%d%d", &i, &j, &next, &hasv)
       || i != d || j != n || next != ab->model) {
      fprintf(stderr, "first line is corrupted:\n%s", s);
      goto ERR;
    }
  }

  if (sizeof(double) == sizeof(real))
    fmt = "%lf%n";
  else
    fmt = "%f%n";
  for (i = 0; i < n; i++) {
    fgets(s, sizeof s, fp);
    if (strlen(s) < 10) goto ERR;
    for (p = s, j = 0; j < d; j++, p += next) {
      if (1 != sscanf(p, fmt, x+i*d+j, &next)) {
        fprintf(stderr, "Error reading i = %d, j = %d\n", i, j);
        goto ERR;
      }
    }
    if (hasv) {
      vi = (v != NULL) ? (v + i*d) : vtmp;
      for (j = 0; j < d; j++, p += next) {
        if (1 != sscanf(p, fmt, vi+j, &next)) {
          fprintf(stderr, "reading i = %d, j = %d\n", i, j);
          goto ERR;
        }
      }
    }
    if (1 != sscanf(p, "%d", &j) || j != ab->type[i]) {
      fprintf(stderr, "bad type on i = %d, j = %d\n", i, j);
      goto ERR;
    }
  }
  fclose(fp);
  return 0;

ERR:
  fprintf(stderr, "position file [%s] appears to be broken on line %d!\n%s\n", fname, i, s);
  fclose(fp);
  return 1;
}

/* initialize a straight chain */
int ab_initpos(abpro_t *ab)
{
  int i, j;

  for (i = 0; i < ab->n; i++) {
    ab->x[i*ab->d] = i - (ab->n - 1)*0.5f;
    for (j = 1; j < ab->d; j++) ab->x[i*ab->d + j] = 0;
  }
  return 0;
}

static int ab_shake2d(abpro_t *ab, const real *x0, real *x1, int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  real dr[2], g, r2;

  for (i = 0; i < n-1; i++)
    rv2_diff(ab->dx + i*2, x0 + (i+1)*2, x0 + i*2);

  for (trial = 0; trial < itmax; trial++) {
    for (again = 0, i = 0; i < n-1; i++) {
      r2 = rv2_sqr(rv2_diff(dr, x1 + (i+1)*2, x1 + i*2));
      if (r2 > 10.0) return 1; /* too large, impossible to correct */

      if (fabs(r2-1) > tol) {
        again = 1;

        g = rv2_dot(dr, ab->dx + i*2);
        if (fabs(g) < 0.1) return 2; /* inner product too small */
        g = (1-r2)/(4*g);
        rv2_sinc(x1 + i*2,     ab->dx + i*2, -g);
        rv2_sinc(x1 + (i+1)*2, ab->dx + i*2,  g);
      }
    }
    if (!again) break;
  }

  if (trial == itmax) return 3;
  return 0;
}

static int ab_shake3d(abpro_t *ab, const real *x0, real *x1, int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  real dr[3], g, r2;

  for (i = 0; i < n-1; i++)
    rv3_diff(ab->dx + i*3, x0 + (i+1)*3, x0 + i*3);

  for (trial = 0; trial < itmax; trial++) {
    for (again = 0, i = 0; i < n-1; i++) {
      r2 = rv3_sqr(rv3_diff(dr, x1 + (i+1)*3, x1 + i*3));
      if (r2 > 10.0) return 1; /* too large, impossible to correct */

      if (fabs(r2-1) > tol) {
        again = 1;

        g = rv3_dot(dr, ab->dx + i*3);
        if (fabs(g) < 0.1) return 2; /* inner product too small */
        g = (1-r2)/(4*g);
        rv3_sinc(x1 + i*3,     ab->dx + i*3, -g);
        rv3_sinc(x1 + (i+1)*3, ab->dx + i*3,  g);
      }
    }
    if (!again) break;
  }
  if (trial == itmax) return 3;
  return 0;
}

/* shake x1 according to x0 */
int ab_shake(abpro_t *ab, const real *x0, real *x1, int itmax, double tol)
{
  if (itmax <= 0) itmax = 10000;
  if (tol <= 0.) tol = 1e-6;
  return (ab->d == 3) ? ab_shake3d(ab, x0, x1, itmax, tol) :
    ab_shake2d(ab, x0, x1, itmax, tol);
}

static int ab_rattle2d(abpro_t *ab, const real *x0, real *v, int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  real dv[2], g;

  for (i = 0; i < n-1; i++)
    rv2_diff(ab->dx + i*2, x0 + (i+1)*2, x0 + i*2);

  for (trial = 0; trial < itmax; trial++) {
    for (again = 0, i = 0; i < n-1; i++) {
      rv2_diff(dv, v + (i+1)*2, v + i*2);
      g = .5f * rv2_dot(ab->dx + i*2, dv);
      if (fabs(g) > tol) {
        again = 1;
        rv2_sinc(v + i*2,     ab->dx + i*2,  g);
        rv2_sinc(v + (i+1)*2, ab->dx + i*2, -g);
      }
    }
    if (!again) break;
  }
  if (trial == itmax) return 1;
  return 0;
}

static int ab_rattle3d(abpro_t *ab, const real *x0, real *v, int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  real dv[3], g;

  for (i = 0; i < n-1; i++)
    rv3_diff(ab->dx + i*3, x0 + (i+1)*3, x0 + i*3);

  for (trial = 0; trial < itmax; trial++) {
    for (again = 0, i = 0; i < n-1; i++) {
      rv3_diff(dv, v + (i+1)*3, v + i*3);
      g = .5f * rv3_dot(ab->dx + i*3, dv);
      if (fabs(g) > tol) {
        again = 1;
        rv3_sinc(v + i*3,     ab->dx + i*3,  g);
        rv3_sinc(v + (i+1)*3, ab->dx + i*3, -g);
      }
    }
    if (!again) break;
  }
  if (trial == itmax) return 1;
  return 0;
}

/* shake v according to x0 */
int ab_rattle(abpro_t *ab, const real *x0, real *v, int itmax, double tol)
{
  if (itmax <= 0) itmax = 10000;
  if (tol <= 0.) tol = 1e-4;
  return (ab->d == 3) ? ab_rattle3d(ab, x0, v, itmax, tol) :
    ab_rattle2d(ab, x0, v, itmax, tol);
}

static int ab_milcshake2d(abpro_t *ab, const real *x0, real *x1, real *v, real dt,
    int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  static real *dl, *dm, *du, *lam, *sig, *x, dx[2];
  register real y;

  if (dl == NULL) {
    if (x0 == NULL) return 0;
    xnew(dl, n); xnew(dm, n); xnew(du, n);
    xnew(lam, n); xnew(sig, n); xnew(x, n*2);
  } else if (x0 == NULL) {
    free(dl); free(dm); free(du);
    free(lam); free(sig); free(x);
    return 0;
  }

  for (i = 0; i < n-1; i++) {
    rv2_diff(ab->dx1 + i*2, x1 + i*2, x1 + (i+1)*2);
    rv2_diff(ab->dx  + i*2, x0 + i*2, x0 + (i+1)*2);
  }

  dm[0] = 2*rv2_dot(ab->dx1, ab->dx);
  du[0] = -rv2_dot(ab->dx1, ab->dx + 2);
  for (i = 1; i < n-1; i++) {
    dl[i] = -rv2_dot(ab->dx1 + i*2, ab->dx + (i-1)*2);
    dm[i] = 2*rv2_dot(ab->dx1 + i*2, ab->dx + i*2);
    du[i] = -rv2_dot(ab->dx1 + i*2, ab->dx + (i+1)*2); /* wrong for i==N-2, but it doesn't hurt */
  }
  for (i = 0; i < n-1; i++)
    sig[i] = 0.5f*(1 - rv2_sqr(ab->dx1 + i*2));

  /* LU decompose D matrix */
  if (fabs(dm[0]) <= 0.) return 1;
  for (i = 1; i < n-1; i++) {
    dm[i] -= dl[i] * (du[i-1] /= dm[i-1]);
    if (fabs(dm[i]) <= 0.) return i+1;
  }

  for (trial = 1; trial <= itmax; trial++) {
    lam[0] = sig[0]/dm[0];
    for (i = 1; i < n-1; i++) /* solving L v = sig */
      lam[i] = (sig[i]-dl[i]*lam[i-1])/dm[i];
    for (i = n-2; i > 0; i--) /* solving U lam = v */
      lam[i-1] -= du[i-1]*lam[i];

    /* if(tridag(n-1, dl,dm,du, lam, sig) != 0) return 1; */
    memcpy(x, x1, 2*n*sizeof(real));
    /* update the new position */
    for (i = 0; i < n-1; i++) {
      rv2_sinc(x + 2*i,     x0 + i*2,  lam[i]);
      rv2_sinc(x + 2*(i+1), x0 + i*2, -lam[i]);
    }

    /* calcualte the maximal error */
    for (again = 0, i = 0; i < n-1; i++) {
      rv2_diff(dx, x + i*2, x + (i+1)*2);
      y = 1 - rv2_sqr(dx);
      if (fabs(y) > tol) again = 1;
      sig[i] += y*.5f;
    }
    if (!again) break;
  }

  memcpy(x1, x, 2*n*sizeof(real));
  if (v != NULL) { /* correct velocities */
    for (i = 0; i < n-1; i++) {
      rv2_sinc(v + 2*i,     x0 + i*2,  lam[i]/dt);
      rv2_sinc(v + 2*(i+1), x0 + i*2, -lam[i]/dt);
    }
  }

  return (trial == itmax);
}

static int ab_milcshake3d(abpro_t *ab, const real *x0, real *x1, real *v, real dt,
    int itmax, double tol)
{
  int i, again, trial, n = ab->n;
  static real *dl, *dm, *du, *lam, *sig, *x, dx[3];
  register real y;

  if (dl == NULL) {
    if (x0 == NULL) return 0;
    xnew(dl, n); xnew(dm, n); xnew(du, n);
    xnew(lam, n); xnew(sig, n); xnew(x, n*3);
  } else if (x0 == NULL) {
    free(dl); free(dm); free(du);
    free(lam); free(sig); free(x);
    return 0;
  }

  for (i = 0; i < n-1; i++) {
    rv3_diff(ab->dx1 + i*3, x1 + i*3, x1 + (i+1)*3);
    rv3_diff(ab->dx  + i*3, x0 + i*3, x0 + (i+1)*3);
  }

  dm[0] = 2*rv3_dot(ab->dx1, ab->dx);
  du[0] = -rv3_dot(ab->dx1, ab->dx + 3);
  for (i = 1; i < n-1; i++) {
    dl[i] = -rv3_dot(ab->dx1 + i*3, ab->dx + (i-1)*3);
    dm[i] = 2*rv3_dot(ab->dx1 + i*3, ab->dx + i*3);
    du[i] = -rv3_dot(ab->dx1 + i*3, ab->dx + (i+1)*3);
  }
  for (i = 0; i < n-1; i++)
    sig[i] = 0.5f*(1 - rv3_sqr(ab->dx1 + i*3));

  /* LU decompose D matrix */
  if (fabs(dm[0]) <= 0.) return 1;
  for (i = 1; i < n-1; i++) {
    dm[i] -= dl[i] * (du[i-1] /= dm[i-1]);
    if (fabs(dm[i]) <= 0.) return i+1;
  }

  for (trial = 1; trial <= itmax; trial++) {
    lam[0] = sig[0]/dm[0];
    for (i = 1; i < n-1; i++) /* solving L v = sig */
      lam[i] = (sig[i]-dl[i]*lam[i-1])/dm[i];
    for (i = n-2; i > 0; i--) /* solving U lam = v */
      lam[i-1] -= du[i-1]*lam[i];

    memcpy(x, x1, 3*n*sizeof(real));
    /* update the new position */
    for (i = 0; i < n-1; i++) {
      rv3_sinc(x + 3*i,     x0 + i*3,  lam[i]);
      rv3_sinc(x + 3*(i+1), x0 + i*3, -lam[i]);
    }

    /* calcualte the maximal error */
    for (again = 0, i = 0; i < n-1; i++) {
      rv3_diff(dx, x + i*3, x + (i+1)*3);
      y = 1 - rv3_sqr(dx);
      if (fabs(y) > tol) again = 1;
      sig[i] += y*.5f;
    }
    if (!again) break;
  }

  memcpy(x1, x, 3*n*sizeof(real));
  if (v != NULL) { /* correct velocities */
    for (i = 0; i < n-1; i++) {
      rv3_sinc(v + 3*i,     x0 + i*3,  lam[i]/dt);
      rv3_sinc(v + 3*(i+1), x0 + i*3, -lam[i]/dt);
    }
  }

  return (trial == itmax);
}

int ab_milcshake(abpro_t *ab, const real *x0, real *x1, real *v, real dt,
    int itmax, double tol)
{
  return (ab->d == 3) ?
    ab_milcshake3d(ab, x0, x1, v, dt, itmax, tol) :
    ab_milcshake2d(ab, x0, x1, v, dt, itmax, tol);
}

static real ab_energy2dm1(abpro_t *ab, const real *r, int soft)
{
  int i, j, n = ab->n;
  real dx[2], dr, dr2, dr6, U = 0;

  for (i = 0; i < n - 1; i++)
    rv2_diff(ab->dx + 2*i, r + 2*(i+1), r + 2*i);

  for (i = 0; i < n - 2; i++)
    U += 1 - rv2_dot(ab->dx + 2*(i+1), ab->dx + 2*i);
  U *= .25f;

  for (i = 0; i < n - 2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv2_sqr( rv2_diff(dx, r + 2*j, r + 2*i) );
      if (soft && dr2 < 1.f) {
        dr = (real) sqrt(dr2);
        printf("%d %d dr = %g\n", i, j, dr);
        U += (52 - 48*dr) - ab->clj[ab->type[i]][ab->type[j]]*(28 - 24*dr);
      } else {
        dr2 = 1/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4*dr6*(dr6 - ab->clj[ab->type[i]][ab->type[j]]);
      }
    }
  }
  return U;
}

static real ab_energy3dm1(abpro_t *ab, const real *r, int soft)
{
  int i, j, n = ab->n;
  real dx[3], dr, dr2, dr6, U = 0;

  for (i = 0; i < n - 1; i++)
    rv3_diff(ab->dx + 3*i, r + 3*(i+1), r + 3*i);

  for (i = 0; i < n - 2; i++)
    U += 1.f - rv3_dot(ab->dx + 3*(i+1), ab->dx + 3*i);
  U *= 0.25f;

  for (i = 0; i < n - 2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_sqr( rv3_diff(dx, r + 3*j, r + 3*i) );
      if (soft && dr2 < 1.f) {
        dr = (real) sqrt(dr2);
        U += (52 - 48*dr) - ab->clj[ab->type[i]][ab->type[j]]*(28 - 24*dr);
      } else {
        dr2 = 1/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4*dr6*(dr6 - ab->clj[ab->type[i]][ab->type[j]]);
      }
    }
  }
  return U;
}

static real ab_energy3dm2(abpro_t *ab, const real *r, int soft)
{
  int i, j, n = ab->n;
  real dx[3], dr2, dr6, U = 0;

  for (i = 0; i < n - 1; i++)
    rv3_diff(ab->dx + 3*i, r + 3*(i+1), r + 3*i);

  for (i = 1; i < n-1; i++)
    U += rv3_dot(ab->dx + 3*i, ab->dx + 3*(i-1));

  for (i = 1; i < n-2; i++)
    U -= .5f * rv3_dot(ab->dx + 3*(i+1), ab->dx + 3*(i-1));

  for (i = 0; i < n-2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_sqr( rv3_diff(dx, r + 3*j, r + 3*i) );
      if (soft && dr2 < 1.f) {
        U += ab->clj[ab->type[i]][ab->type[j]]*(ab->sla - ab->slb*(real)sqrt(dr2));
      } else {
        dr2 = 1.f/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4.f*ab->clj[ab->type[i]][ab->type[j]]*dr6*(dr6 - 1.f);
      }
    }
  }
  return U;
}

real ab_energy(abpro_t *ab, const real *r, int soft)
{
  if (ab->model == 2)
    return ab_energy3dm2(ab, r, soft);
  else if (ab->d == 3)
    return ab_energy3dm1(ab, r, soft);
  else
    return ab_energy2dm1(ab, r, soft);
}

static real ab_force2dm1(abpro_t *ab, real *f, const real *r, int soft)
{
  int i, j, n = ab->n;
  real c, ff, dr2, dr6, U = 0.f;
  real *dxp, *dxm, dx[2];

  for (i = 0; i < n; i++) rv2_zero(f + 2*i);

  for (i = 0; i < n - 1; i++)
    rv2_diff(ab->dx + 2*i, r + 2*(i+1), r + 2*i);

  for (i = 1; i < n-1; i++) {
    dxm = (dxp = ab->dx + 2*i) - 2;
    U += 1.f - rv2_dot(dxp, dxm);
    rv2_sinc(f + 2*(i-1), dxp, -.25f);
    rv2_sinc(f + 2*i,     dxp,  .25f);
    rv2_sinc(f + 2*i,     dxp, -.25f);
    rv2_sinc(f + 2*(i+1), dxp,  .25f);
  }
  U *= 0.25f;

  for (i = 0; i < n-2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv2_sqr( rv2_diff(dx, r + 2*j, r + 2*i) );
      c = ab->clj[ab->type[i]][ab->type[j]];

      if (soft && dr2 < 1.) {
        dr2 = (real) sqrt(dr2);
        U += (52.f - 28.f*c) - 24.f*dr2*(2.f-c);
        ff = 24.f*(2.f - c)/dr2;
      } else {
        dr2 = 1.f/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4.f*dr6*(dr6 - c);
        ff = 24.f*dr2*dr6*(dr6*2.f - c);
      }
      rv2_sinc(f + 2*i, dx, -ff);
      rv2_sinc(f + 2*j, dx, +ff);
    }
  }
  return U;
}

static real ab_force3dm1(abpro_t *ab, real *f, const real *r, int soft)
{
  int i, j, n = ab->n;
  real c, ff, dr2, dr6, U = 0.f;
  real *dxp, *dxm, dx[3];

  for (i = 0; i < n; i++) rv3_zero(f + 3*i);

  for (i = 0; i < n - 1; i++)
    rv3_diff(ab->dx + 3*i, r + 3*(i+1), r + 3*i);

  for (i = 1; i < n-1; i++) {
    dxm = (dxp = ab->dx + 3*i) - 3;
    U += 1.f - rv3_dot(dxp, dxm);
    rv3_sinc(f + 3*(i-1), dxp, -.25f);
    rv3_sinc(f + 3*i,     dxp,  .25f);
    rv3_sinc(f + 3*i,     dxp, -.25f);
    rv3_sinc(f + 3*(i+1), dxp,  .25f);
  }
  U *= 0.25f;

  for (i = 0; i < n-2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_sqr( rv3_diff(dx, r + 3*j, r + 3*i) );
      c = ab->clj[ab->type[i]][ab->type[j]];

      if (soft && dr2 < 1.) {
        dr2 = (real) sqrt(dr2);
        U += (52.f - 28.f*c) - 24.f*dr2*(2.f-c);
        ff = 24.f*(2.f - c)/dr2;
      } else {
        dr2 = 1.f/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4.f*dr6*(dr6 - c);
        ff = 24.f*dr2*dr6*(dr6*2.f - c);
      }
      rv3_sinc(f + 3*i, dx, -ff);
      rv3_sinc(f + 3*j, dx, +ff);
    }
  }
  return U;
}


static real ab_force3dm2(abpro_t *ab, real *f, const real *r, int soft)
{
  real ff, dr2, dr6, U = 0;
  real *dxm, *dxp, dx[3], c;
  int i, j, n = ab->n;

  for (i = 0; i < n; i++) rv3_zero(f + 3*i);

  for (i = 0; i < n - 1; i++)
    rv3_diff(ab->dx + 3*i, r + 3*(i+1), r + 3*i);

  for (i = 1; i < n-1; i++) {
    dxm = (dxp = ab->dx + 3*i) - 3;
    rv3_inc(f + (i-1)*3, dxp);
    rv3_dec(f + i*3,     dxp);
    rv3_inc(f + i*3,     dxm);
    rv3_dec(f + (i+1)*3, dxm);
    U += rv3_dot(dxp, dxm);
  }

  for (i = 1; i < n-2; i++) {
    dxp = ab->dx + 3*(i+1);
    dxm = ab->dx + 3*(i-1);
    rv3_sinc(f + 3*(i-1), dxp, -.5f);
    rv3_sinc(f + 3*i,     dxp,  .5f);
    rv3_sinc(f + 3*(i+1), dxp, -.5f);
    rv3_sinc(f + 3*(i+2), dxp,  .5f);
    U -= .5f*rv3_dot(dxp, dxm);
  }

  for (i = 0; i < n - 2; i++) {
    for (j = i+2; j < n; j++) {
      c = ab->clj[ab->type[i]][ab->type[j]];
      dr2 = rv3_sqr( rv3_diff(dx, r + 3*j, r + 3*i) );
      if (soft && dr2 < 1.f) {
        dr2 = (real) sqrt(dr2);
        U += c*(ab->sla - ab->slb*dr2);
        ff = (ab->slb*c)/dr2;
      } else {
        dr2 = 1.f/dr2;
        dr6 = dr2*dr2*dr2;
        U += 4.f*c*dr6*(dr6 - 1.f);
        ff = 48.f*c*dr2*dr6*(dr6 - .5f);
      }
      rv3_sinc(f + 3*i, dx, -ff);
      rv3_sinc(f + 3*j, dx,  ff);
    }
  }

  return U;
}

/* compute force f */
real ab_force(abpro_t *ab, real *f, real *r, int soft)
{
  if (ab->model == 2)
    return ab_force3dm2(ab, f, r, soft);
  else if (ab->d == 3)
    return ab_force3dm1(ab, f, r, soft);
  else
    return ab_force2dm1(ab, f, r, soft);
}

/* minimizes the energy of a given configuration.
   The minimized configuration is saved in ab->lmx
   When a lowest energy configuration is found, the result is
   saved to global variable ab->xmin, with ab->emin updated. */
real ab_localmin(abpro_t *ab, const real *r, int itmax, double tol)
{
  int t, i, j, id, n, d;
  real up, u0, u = 0, step = 0.02, del, mem = 1;
  static real *x[2], *f[2], *v;
  static double ncnt = 0, cnt = 0;
  const real DELMAX = 0.20f;

  if (ab == NULL && v != NULL) {
    free(x[0]); free(x[1]); free(f[0]); free(f[1]); free(v);
    return 0;
  }
  if (itmax <= 0) itmax = 10000;
  if (tol <= 0.) tol = 1e-12;
  n = ab->n;
  d = ab->d;
  if (v == NULL) {
    xnew(x[0], n*d*sizeof(real));
    xnew(x[1], n*d*sizeof(real));
    xnew(f[0], n*d*sizeof(real));
    xnew(f[1], n*d*sizeof(real));
    xnew(v, n*d*sizeof(real));
  }
  /* to make a working copy */
  memcpy(x[id = 0], r, n*d*sizeof(real));
  up = ab_force(ab, f[id], x[id], 0);
  memset(v, 0, n*d*sizeof(real));
  u0 = up;

  for (t = 1; t <= itmax; t++) {
    for (i = 0; i < n; i++)
      for (j = 0; j < d; j++) {
        del = v[i*d+j] = v[i*d+j]*mem+f[id][i*d+j]*step;
        if (del > DELMAX) del = DELMAX; else if (del < -DELMAX) del = -DELMAX;
        x[!id][i*d+j] = x[id][i*d+j]+del;
      }

    if (ab_shake(ab, x[id], x[!id], 0, 0.) != 0) goto SHRINK;
    u = ab_force(ab, f[!id], x[!id], 0);
    if (u > up) { mem = 0; goto SHRINK; }

    id = !id;
    if (up - u < tol) break;
    up = u;
    mem = 0.9;
    step *= 1.1;
    continue;

SHRINK:
    step *= 0.5;
  }
  if (t > itmax) fprintf(stderr, "local minimum failed to converge.\n");

  cnt++, ncnt += t;
  if (fmod(cnt, 100) < 0.1)
    fprintf(stderr, "abmin: cnt=%g aver=%g u: %g => %g\n", cnt, 1.0*ncnt/cnt, u0, u);

  memcpy(ab->lmx, x[id], n*d*sizeof(real));
  if (u < ab->emin) {
    ab->emin = u;
    memcpy(ab->xmin, x[id], n*d*sizeof(real));
    ab_writepos(ab, ab->xmin, NULL, "abmin.pos");
  }

  return u;
}

#endif

