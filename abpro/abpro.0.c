#include "util.c"
#include "rng.c"

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
abpro_t *ab_open(int seqid, int d, int model, real randdev)
{
  abpro_t *ab;
  int i, nd;
  double x;

  die_if (d == 2 && model != 1, "%dd only for model 1", d);
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

  ab->seqid = seqid;
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

  /* number of degrees of freedom */
  ab->dof = (ab->d == 2) ? (ab->n - 2) : (2*ab->n - 5);
  //printf("n = %3d, dof = %3d: ", ab->n, ab->dof);
  //for (i = 0; i < ab->n; i++) printf("%c", ab->type[i]+'A'); printf("\n");

  nd = ab->n * ab->d;
  xnew(ab->x, nd);
  xnew(ab->x1, nd);
  xnew(ab->dx, nd);
  xnew(ab->v, nd);
  xnew(ab->f, nd);
  xnew(ab->lmx, nd);
  xnew(ab->xmin, nd);

  ab_initpos(ab, ab->x, randdev);
  ab->emin = ab->epot = ab_force(ab, ab->f, ab->x, 0);
  return ab;
}

/* initialize an almost straight chain,
 * randomness given by del */
int ab_initpos(abpro_t *ab, real *x, real del)
{
  int i, j;
  real dx[3];

  for (j = 0; j < ab->d; j++) ab->x[j] = 0;
  for (i = 0; i < ab->n - 1; i++) {
    for (j = 0; j < ab->d; j++) {
      dx[j] = (2.f*rand()/RAND_MAX - 1)*del + ((j == 0) ? 1.f : 0.f);
    }
    if (ab->d == 3) {
      rv3_normalize(dx);
      rv3_add(x + (i+1)*ab->d, x + i*ab->d, dx);
    } else {
      rv2_normalize(dx);
      rv2_add(x + (i+1)*ab->d, x + i*ab->d, dx);
    }
  }
  die_if (ab_checkconn(ab, x, 0) != 0, "initpos failed, with del = %g\n", del);
  ab_shiftcom(ab, x);
  return 0;
}

/* close ab */
void ab_close(abpro_t *ab)
{
  if (ab) {
    ab_milcshake(ab, NULL, NULL, NULL, 0, 0, 0, 0);
    ab_milcrattle(ab, NULL, NULL);
    ab_localmin(ab, NULL, 0, 0., 0, 0., 0);
    free(ab->type);
    free(ab->x);
    free(ab->x1);
    free(ab->dx);
    free(ab->v);
    free(ab->f);
    free(ab->lmx);
    free(ab->xmin);
    free(ab);
  }
}

/* check connectivity */
int ab_checkconn(abpro_t *ab, const real *x, double tol)
{
  int i, d = ab->d;
  real r;

  if (tol <= 0.) tol = 1e-3;
  for (i = 0; i < ab->n-1; i++) {
    if (d == 3) {
      r = rv3_dist(x + i*3, x + (i+1)*3);
    } else {
      r = rv2_dist(x + i*2, x + (i+1)*2);
    }
    if (fabs(r-1) > tol) {
      fprintf(stderr, "link (%d,%d) is broken, r = %g\n", i, i+1, r);
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

/* annihilate angular momentum 2d */
static void ab_rmangular2d(abpro_t *ab, crv2_t *x, rv2_t *v)
{
  int i, n = ab->n;
  real am, r2, rp[2];

  for (am = r2 = 0.f, i = 0; i < n; i++) {
    am += rv2_cross(x[i], v[i]);
    r2 += rv2_sqr(x[i]);
  }
  am = -am/r2;
  for (i = 0; i < n; i++) {
    rp[0] = -am*x[i][1];
    rp[1] =  am*x[i][0];
    rv2_inc(v[i], rp);
  }
}

/* annihilate angular momentum 3d
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  / 
 * use
 *    v = c X r
 *   */
static void ab_rmangular3d(abpro_t *ab, crv3_t *x, rv3_t *v)
{
  int i, n = ab->n;
  real ang[3], am[3], dv[3], mat[3][3], inv[3][3];
  const real *xi;
  real xx = 0.f, yy = 0.f, zz = 0.f, xy = 0.f, zx = 0.f, yz = 0.f;

  rv3_zero(am);
  for (i = 0; i < n; i++) {
    rv3_cross(ang, x[i], v[i]);
    rv3_inc(am, ang);
    xi = x[i];
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  mat3_inv(inv, mat);
  ang[0] = -rv3_dot(inv[0], am);
  ang[1] = -rv3_dot(inv[1], am);
  ang[2] = -rv3_dot(inv[2], am);
  /* ang is the solution of M^(-1) * I */
  for (i = 0; i < n; i++) {
    rv3_cross(dv, ang, x[i]);
    rv3_inc(v[i], dv);
  }
}

/* shift center of x to the origin, 
 * remove center velocity and angular momentum */
void ab_rmcom(abpro_t *ab, real *x, real *v)
{
  ab_shiftcom(ab, x);
  ab_shiftcom(ab, v);
  /* remove angular momentum */
  if (ab->d == 2) ab_rmangular2d(ab, (crv2_t *)x, (rv2_t *)v); 
  else ab_rmangular3d(ab, (crv3_t *)x, (rv3_t *)v);
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

  fprintf(fp, "# %d %d %d %d %d\n", d, ab->model, ab->seqid, ab->n, (v != NULL));
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
  int i, j, seq, hasv = 0, next, d = ab->d, n = ab->n;
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
    if (5 != sscanf(s+1, "%d%d%d%d%d", &i, &j, &seq, &next, &hasv)
       || i != d || j != ab->model || seq != ab->seqid || next != n) {
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

static int ab_shake3d(abpro_t *ab, crv3_t *x0, rv3_t *x1, rv3_t *v, real dt,
    int itmax, double tol, int verbose)
{
  int i, again, it, n = ab->n;
  real dxi[3], g, r2, r2bad;
  rv3_t *dx0 = (rv3_t *)ab->dx;
  const real glow = .5, r2max = 4.0;

  for (i = 0; i < n-1; i++)
    rv3_diff(dx0[i], x0[i+1], x0[i]);

  for (it = 0; it < itmax; it++) {
    for (again = 0, i = 0; i < n-1; i++) {
      r2 = rv3_sqr(rv3_diff(dxi, x1[i+1], x1[i]));
      if (r2 > r2max) { /* too large, impossible to correct */
        if (verbose)
          fprintf(stderr, "shake: large distance %d-%d, %g\n", i,i+1, sqrt(r2));
        r2 = r2max; 
      }

      if (fabs(r2-1) > tol) {
        if (!again) { again = 1; r2bad = r2; }

        g = rv3_dot(dxi, dx0[i]);
        if (fabs(g) < glow) { /* inner product too small */
          if (verbose)
            fprintf(stderr, "shake: bad alignment %d-%d, %g, dot = %g\n", i,i+1, sqrt(r2), g);
          g = (g > 0) ? glow : -glow;
        }
        g = (1-r2)/(4*g);
        rv3_sinc(x1[i],   dx0[i], -g);
        rv3_sinc(x1[i+1], dx0[i],  g);
        if (v) { /* add a force of dx/dt */
          rv3_sinc(v[i],   dx0[i], -g/dt);
          rv3_sinc(v[i+1], dx0[i],  g/dt);
        }
      }
    }
    if (!again) break;
  }
  if (it >= itmax) {
    if (verbose) {
      const char *fnf = "shakefail.pos";
      fprintf(stderr, "shake: failed after %d iter. r = 1%+g, see %s\n", it, sqrt(r2bad)-1, fnf);
      ab_writepos(ab, (real *)x1, NULL, fnf);
    }
    return -1;
  }
  return 0;
}

/* shake x1 according to x0 */
int ab_shake(abpro_t *ab, const real *x0, real *x1, real *v, real dt,
    int itmax, double tol, int verbose)
{
  if (itmax <= 0) itmax = 3000;
  if (tol <= 0.) tol = (sizeof(real) == sizeof(double)) ? 1e-6 : 1e-4;
  return (ab->d == 3) ? 
    ab_shake3d(ab, (crv3_t *)x0, (rv3_t *)x1, (rv3_t *)v, dt, itmax, tol, verbose) :
    ab_shake2d(ab, (crv2_t *)x0, (rv2_t *)x1, (rv2_t *)v, dt, itmax, tol, verbose);
}

static int ab_rattle3d(abpro_t *ab, crv3_t *x0, rv3_t *v, 
    int itmax, double tol, int verbose)
{
  int i, again, it, n = ab->n;
  real dv[3], g, rvbad;
  rv3_t *dx = (rv3_t *)ab->dx;

  for (i = 0; i < n-1; i++)
    rv3_diff(dx[i], x0[i+1], x0[i]);

  for (it = 0; it < itmax; it++) {
    for (again = 0, i = 0; i < n-1; i++) {
      rv3_diff(dv, v[i+1], v[i]);
      g = .5f * rv3_dot(dx[i], dv);
      if (fabs(g) > tol) {
        if (!again) { again = 1; rvbad = g; }
        rv3_sinc(v[i],   dx[i],  g);
        rv3_sinc(v[i+1], dx[i], -g);
      }
    }
    if (!again) break;
  }
  if (it >= itmax) {
    if (verbose) {
      const char *fnf = "rattlefail.pos";
      fprintf(stderr, "rattle: failed after %d iter. rv = %+g, see %s\n", it, rvbad, fnf);
      ab_writepos(ab, (real *)x0, (real *)v, fnf);
    }    
    return -1;
  }
  return 0;
}

/* shake v according to x0 */
int ab_rattle(abpro_t *ab, const real *x0, real *v, int itmax, double tol, int verbose)
{
  if (itmax <= 0) itmax = 3000;
  if (tol <= 0.) tol = 1e-4;
  return (ab->d == 3) ? 
    ab_rattle3d(ab, (crv3_t *)x0, (rv3_t *)v, itmax, tol, verbose) :
    ab_rattle2d(ab, (crv2_t *)x0, (rv2_t *)v, itmax, tol, verbose);
}

static int ab_milcshake3d(abpro_t *ab, crv3_t *x0, rv3_t *x1, rv3_t *v, real dt,
    int itmax, double tol, int verbose)
{
  int i, again, it, n = ab->n, nl;
  static real *dl, *dm, *du, *lam, *rhs;
  static rv3_t *x, *dx0, *dx1;
  real y;

  if (dl == NULL) {
    if (x0 == NULL) return 0;
    xnew(dl, n); xnew(dm, n); xnew(du, n);
    xnew(lam, n); xnew(rhs, n); 
    xnew(dx0, n); xnew(dx1, n); xnew(x, n);
  } else if (x0 == NULL) {
    free(dl); free(dm); free(du);
    free(lam); free(rhs);
    free(dx0); free(dx1); free(x);
    return 0;
  }

  nl = n - 1;
  for (i = 0; i < nl; i++) {
    rv3_diff(dx1[i], x1[i], x1[i+1]);
  }
  for (i = 0; i < nl; i++) {
    rv3_diff(dx0[i], x0[i], x0[i+1]);
  }

  /* dm[0..nl-1], du[0..nl-2], dl[1..nl-1] */
  dm[0] =  4*rv3_dot(dx0[0], dx1[0]);
  du[0] = -2*rv3_dot(dx0[1], dx1[0]);
  for (i = 1; i < nl; i++) {
    dl[i] = -2*rv3_dot(dx1[i], dx0[i-1]);
    dm[i] =  4*rv3_dot(dx1[i], dx0[i]);
    du[i] = -2*rv3_dot(dx1[i], dx0[i+1]); /* no dx0[nl], but doesn't matter */ 
  }
  for (i = 0; i < nl; i++)
    rhs[i] = 1 - rv3_sqr(dx1[i]);

  /* solve matrix equation D lam = rhs
   * first LU decompose D; 
   * U --> du with diagonal being unity; 
   * L --> dm and dl with dl unchanged */
  if (fabs(dm[0]) < 1e-6) return 1;
  for (i = 1; i < nl; i++) {
    dm[i] -= dl[i] * (du[i-1] /= dm[i-1]);
    if (fabs(dm[i]) < 1e-6) return i+1;
  }

  for (it = 1; it <= itmax; it++) {
    lam[0] = rhs[0]/dm[0];
    for (i = 1; i < nl; i++) /* solving L v = rhs */
      lam[i] = (rhs[i] - dl[i]*lam[i-1])/dm[i];
    for (i = nl - 1; i > 0; i--) /* solving U lam = v */
      lam[i-1] -= du[i-1]*lam[i];

    rv3_ncopy(x, x1, n);
    /* update the new position */
    for (i = 0; i < nl; i++) {
      rv3_sinc(x[i],   dx0[i],  lam[i]);
      rv3_sinc(x[i+1], dx0[i], -lam[i]);
    }

    /* calcualte the maximal error */
    for (again = 0, i = 0; i < nl; i++) {
      y = 1 - rv3_dist2(x[i], x[i+1]);
      if (fabs(y) > tol) again = 1;
      rhs[i] += y;
    }
    if (!again) break;
  }

  rv3_ncopy(x1, x, n);
  if (v != NULL) { /* correct velocities */
    for (i = 0; i < n-1; i++) {
      rv3_sinc(v[i],   dx0[i],  lam[i]/dt);
      rv3_sinc(v[i+1], dx0[i], -lam[i]/dt);
    }
  }

  if (it >= itmax) {
    if (verbose) {
      const char *fnf = "shakefail.pos";
      fprintf(stderr, "milcshake: failed after %d iter. see %s\n", it, fnf);
      ab_writepos(ab, (real *)x1, NULL, fnf);
    }    
    return -1;
  }
  return 0;
}

/* MILC shake, make |dr| = 1
 * for a random config., about 30~40% faster than shake
 * but slower than shake for near-minimum config.  */
int ab_milcshake(abpro_t *ab, const real *x0, real *x1, real *v, real dt,
    int itmax, double tol, int verbose)
{
  if (itmax <= 0) itmax = 3000;
  if (tol <= 0.) tol = (sizeof(real) == sizeof(double)) ? 1e-6 : 1e-4;
  return (ab->d == 3) ?
    ab_milcshake3d(ab, (crv3_t *)x0, (rv3_t *)x1, (rv3_t *)v, dt, itmax, tol, verbose) :
    ab_milcshake2d(ab, (crv2_t *)x0, (rv2_t *)x1, (rv2_t *)v, dt, itmax, tol, verbose);
}

static int ab_milcrattle3d(abpro_t *ab, crv3_t *x, rv3_t *v)
{
  int i, n = ab->n, nl;
  static real *dl, *dm, *du, *lam, *rhs;
  static rv3_t *dx, *dv;

  if (dl == NULL) {
    if (x == NULL) return 0;
    xnew(dl, n); xnew(dm, n); xnew(du, n);
    xnew(lam, n); xnew(rhs, n); 
    xnew(dx, n); xnew(dv, n);
  } else if (x == NULL) {
    free(dl); free(dm); free(du);
    free(lam); free(rhs);
    free(dx); free(dv);
    return 0;
  }

  nl = n - 1;
  for (i = 0; i < nl; i++) {
    rv3_diff(dx[i], x[i], x[i+1]);
    rv3_diff(dv[i], v[i], v[i+1]);
  }

  /* dm[0..nl-1], du[0..nl-2], dl[1..nl-1] */
  dm[0] = 2*rv3_dot(dx[0], dx[0]);
  du[0] =  -rv3_dot(dx[1], dx[0]);
  for (i = 1; i < nl; i++) {
    dl[i] =  -rv3_dot(dx[i], dx[i-1]);
    dm[i] = 2*rv3_dot(dx[i], dx[i]);
    du[i] =  -rv3_dot(dx[i], dx[i+1]); /* no dx[nl], but doesn't matter */ 
  }
  for (i = 0; i < nl; i++)
    rhs[i] = -rv3_dot(dv[i], dx[i]);

  /* solve matrix equation D lam = rhs
   * first LU decompose D; 
   * U --> du with diagonal being unity; 
   * L --> dm and dl with dl unchanged */
  if (fabs(dm[0]) < 1e-6) return 1;
  for (i = 1; i < nl; i++) {
    dm[i] -= dl[i] * (du[i-1] /= dm[i-1]);
    if (fabs(dm[i]) < 1e-6) return i+1;
  }

  lam[0] = rhs[0]/dm[0];
  for (i = 1; i < nl; i++) /* solving L v = rhs */
    lam[i] = (rhs[i] - dl[i]*lam[i-1])/dm[i];
  for (i = nl - 1; i > 0; i--) /* solving U lam = v */
    lam[i-1] -= du[i-1]*lam[i];

  /* update the new position */
  for (i = 0; i < nl; i++) {
    rv3_sinc(v[i],   dx[i],  lam[i]);
    rv3_sinc(v[i+1], dx[i], -lam[i]);
  }
  return 0;
}

/* MILC rattle, make dr.v = 0 */
int ab_milcrattle(abpro_t *ab, const real *x, real *v) 
{
  return (ab->d == 3) ?
    ab_milcrattle3d(ab, (crv3_t *)x, (rv3_t *)v) :
    ab_milcrattle2d(ab, (crv2_t *)x, (rv2_t *)v);
}

static real ab_energy3dm1(abpro_t *ab, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  real dr, dr2, dr6, U = 0;
  rv3_t *dx = (rv3_t *)ab->dx;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 0; i < n - 2; i++)
    U += 1.f - rv3_dot(dx[i+1], dx[i]);
  U *= 0.25f;

  for (i = 0; i < n - 2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_dist2(r[j], r[i]);
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

static real ab_energy3dm2(abpro_t *ab, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  real dr2, dr6, U = 0;
  rv3_t *dx = (rv3_t *)ab->dx;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 1; i < n-1; i++)
    U += rv3_dot(dx[i], dx[i-1]);

  for (i = 1; i < n-2; i++)
    U -= .5f * rv3_dot(dx[i+1], dx[i-1]);

  for (i = 0; i < n-2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_dist2(r[j], r[i]);
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
    return ab_energy3dm2(ab, (crv3_t *)r, soft);
  else if (ab->d == 3)
    return ab_energy3dm1(ab, (crv3_t *)r, soft);
  else
    return ab_energy2dm1(ab, (crv2_t *)r, soft);
}

static real ab_force3dm1(abpro_t *ab, rv3_t *f, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  real c, ff, dr2, dr6, U = 0.f;
  real *dxp, *dxm;
  rv3_t *dx = (rv3_t *)ab->dx, dxi;

  for (i = 0; i < n; i++) rv3_zero(f[i]);

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 1; i < n-1; i++) {
    dxp = dx[i];
    dxm = dx[i-1];
    U += 1.f - rv3_dot(dxp, dxm);
    rv3_sinc(f[i-1], dxp, -.25f);
    rv3_sinc(f[i],   dxp,  .25f);
    rv3_sinc(f[i],   dxm, -.25f);
    rv3_sinc(f[i+1], dxm,  .25f);
  }
  U *= 0.25f;

  for (i = 0; i < n-2; i++) {
    for (j = i+2; j < n; j++) {
      dr2 = rv3_sqr( rv3_diff(dxi, r[j], r[i]) );
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
      rv3_sinc(f[i], dxi, -ff);
      rv3_sinc(f[j], dxi, +ff);
    }
  }
  return U;
}

static real ab_force3dm2(abpro_t *ab, rv3_t *f, crv3_t *r, int soft)
{
  real ff, dr2, dr6, U = 0;
  real *dxm, *dxp, c;
  int i, j, n = ab->n;
  rv3_t *dx = (rv3_t *)ab->dx, dxi;

  for (i = 0; i < n; i++) rv3_zero(f[i]);

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 1; i < n-1; i++) {
    dxp = dx[i];
    dxm = dx[i-1];
    rv3_inc(f[i-1], dxp);
    rv3_dec(f[i],   dxp);
    rv3_inc(f[i],   dxm);
    rv3_dec(f[i+1], dxm);
    U += rv3_dot(dxp, dxm);
  }

  for (i = 1; i < n-2; i++) {
    dxp = dx[i+1];
    dxm = dx[i-1];
    rv3_sinc(f[i-1], dxp, -.5f);
    rv3_sinc(f[i],   dxp,  .5f);
    rv3_sinc(f[i+1], dxm, -.5f);
    rv3_sinc(f[i+2], dxm,  .5f);
    U -= .5f*rv3_dot(dxp, dxm);
  }

  for (i = 0; i < n - 2; i++) {
    for (j = i+2; j < n; j++) {
      c = ab->clj[ab->type[i]][ab->type[j]];
      dr2 = rv3_sqr( rv3_diff(dxi, r[j], r[i]) );
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
      rv3_sinc(f[i], dxi, -ff);
      rv3_sinc(f[j], dxi,  ff);
    }
  }

  return U;
}

/* compute force f */
real ab_force(abpro_t *ab, real *f, const real *r, int soft)
{
  if (ab->d == 2)
    return ab_force2dm1(ab, (rv2_t *)f, (crv2_t *)r, soft);
  else if (ab->model == 1)
    return ab_force3dm1(ab, (rv3_t *)f, (crv3_t *)r, soft);
  else
    return ab_force3dm2(ab, (rv3_t *)f, (crv3_t *)r, soft);
}

/* minimizes the energy of a given configuration.
   The minimized configuration is saved in ab->lmx
   When a lowest energy configuration is found, the result is
   saved to global variable ab->xmin, with ab->emin updated. */
real ab_localmin(abpro_t *ab, const real *r, int itmax, double tol,
    int sh_itmax, double sh_tol, unsigned flags)
{
  int t, i, j, id, n = ab->n, d = ab->d;
  real up, u = 0, step = 0.02, del, mem = 1;
  static real *x[2], *f[2], *v;
  const real DELMAX = 0.20f;

  if (v == NULL) {
    if (r == NULL) return 0;
    xnew(x[0], n*d*sizeof(real));
    xnew(x[1], n*d*sizeof(real));
    xnew(f[0], n*d*sizeof(real));
    xnew(f[1], n*d*sizeof(real));
    xnew(v, n*d*sizeof(real));
  } else if (r == NULL) {
    free(x[0]); free(x[1]); free(f[0]); free(f[1]); free(v);
    return 0;
  }
  if (itmax <= 0) itmax = 10000;
  if (tol <= 0.) tol = 1e-12;
  /* to make a working copy */
  memcpy(x[id = 0], r, n*d*sizeof(real));
  up = ab_force(ab, f[id], x[id], 0);
  memset(v, 0, n*d*sizeof(real));

  for (t = 1; t <= itmax; t++) {
    for (i = 0; i < n; i++)
      for (j = 0; j < d; j++) {
        del = v[i*d+j] = v[i*d+j]*mem+f[id][i*d+j]*step;
        if (del > DELMAX) del = DELMAX; else if (del < -DELMAX) del = -DELMAX;
        x[!id][i*d+j] = x[id][i*d+j]+del;
      }

    if (flags & AB_MILCSHAKE) {
      if (ab_milcshake(ab, x[id], x[!id], NULL, 0., sh_itmax, sh_tol, 0) != 0) goto SHRINK;
    } else {
      if (ab_shake(ab, x[id], x[!id], NULL, 0., sh_itmax, sh_tol, 0) != 0) goto SHRINK;
    }
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
  if (t > itmax && (flags & AB_VERBOSE)) 
    fprintf(stderr, "localmin failed to converge, t = %d.\n", t);

  memcpy(ab->lmx, x[id], n*d*sizeof(real));
  if (u < ab->emin && (flags & AB_LMREGISTER)) {
    ab->emin = u;
    memcpy(ab->xmin, x[id], n*d*sizeof(real));
    if (flags & AB_LMWRITE)
      ab_writepos(ab, ab->xmin, NULL, "abmin.pos");
  }

  return u;
}

/* velocity rescaling */
void ab_vrescale(abpro_t *ab, real tp, real dt)
{
  int i;
  real ekav = .5f*tp*ab->dof, ek1 = ab->ekin, ek2, s;
  double amp;

  amp = 2*sqrt(ek1*ekav*dt/ab->dof);
  ek2 = ek1 + (ekav - ek1)*dt + (real)(amp*grand0());
  s = (real)sqrt(ek2/ek1);
  for (ab->ekin = 0.f, i = 0; i < ab->n*ab->d; i++)
    ab->v[i] *= s;
  ab->ekin = ek2;
  ab->tkin *= s*s;
}

/* kinetic energy */
real ab_ekin(abpro_t *ab)
{
  int i, nd = ab->n * ab->d;

  for (ab->ekin = 0, i = 0; i < nd; i++) 
    ab->ekin += ab->v[i]*ab->v[i];
  ab->tkin = ab->ekin/ab->dof;
  return (ab->ekin *= .5f);
}

static int ab_vv3d(abpro_t *ab, real fscal, real dt, int soft, int milc)
{
  int i, verbose = 1, n = ab->n;
  real dth = .5f*dt*fscal;
  rv3_t *v = (rv3_t *)ab->v, *x = (rv3_t *)ab->x, *x1 = (rv3_t *)ab->x1, *f = (rv3_t *)ab->f;

  for (i = 0; i < n; i++) { /* vv part 1 */
    rv3_sinc(v[i], f[i], dth);
    rv3_lincomb2(x1[i], x[i], v[i], 1, dt);
  }
  if (milc) {
    i = ab_milcshake(ab, ab->x, ab->x1, ab->v, dt, 0, 0., verbose);
  } else {
    i = ab_shake(ab, ab->x, ab->x1, ab->v, dt, 0, 0., verbose);
  }
  die_if (i != 0, "t=%g: shake failed\n", ab->t);
  rv3_ncopy(x, x1, n);

  ab->epot = ab_force(ab, ab->f, ab->x, soft); /* calculate force */
  
  for (i = 0; i < n; i++) { /* vv part 2 */
    rv3_sinc(v[i], f[i], dth);
  }
  if (milc) {
    i = ab_milcrattle(ab, ab->x, ab->v);
  } else {
    i = ab_rattle(ab, ab->x, ab->v, 0, 0., verbose);
  }

  ab_ekin(ab);

  die_if (i != 0, "t=%g: failed rattle\n", ab->t);
  ab->t += dt;
  return 0;
}

/* one step of velocity verlet integrator */
int ab_vv(abpro_t *ab, real fscal, real dt, unsigned flags)
{
  int soft = (flags & AB_SOFTFORCE), milc = (flags & AB_MILCSHAKE);
  return (ab->d == 3) ?
    ab_vv3d(ab, fscal, dt, soft, milc) :
    ab_vv2d(ab, fscal, dt, soft, milc);
}

/* Brownian dynamics */
int ab_brownian(abpro_t *ab, real T, real fscal, real dt, unsigned flags)
{
  int soft = (flags & AB_SOFTFORCE), milc = (flags & AB_MILCSHAKE);
  int i, nd = ab->n * ab->d, verbose = 1;
  real amp = (real) sqrt(2*dt*T);

  for (i = 0; i < nd; i++) 
    ab->x1[i] = ab->x[i] + fscal*ab->f[i]*dt + (real)(grand0()*amp);

  if (milc) {
    i = ab_milcshake(ab, ab->x, ab->x1, NULL, 0.f, 0, 0., verbose);
  } else {
    i = ab_shake(ab, ab->x, ab->x1, NULL, 0.f, 0, 0., verbose);
  }
  die_if (i != 0, "t=%g: failed shake\n", ab->t);
  memcpy(ab->x, ab->x1, nd*sizeof(real));

  ab->epot = ab_force(ab, ab->f, ab->x, soft); /* calculate force */
  ab->t += dt;
  return 0;
} 

#endif

