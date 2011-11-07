#ifndef ABPRO_C__
#define ABPRO_C__
/* AB beads protein models */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "abpro.h"

#ifdef _OPENMP
/* compare pair by thread id */
static int ab_prcmp(const void *a, const void *b)
{ return ((abpairid_t *) a)->tid - ((abpairid_t *) b)->tid; }
#endif

/* initialization
 * seqid: 8: 34, 9: 55, 10: 89 */
abpro_t *ab_open(int seqid, int d, int model, real randdev)
{
  abpro_t *ab;
  int i, j, nd;
  double x;
  const int verbose = 0;

  die_if (d == 2 && model != 1, "%dd only for model 1", d);
  die_if (seqid < 0, "bad seqid %d\n", seqid);

  xnew(ab, 1);
  ab->d = d;

  ab->model = model;
  if (model == 1) { /* initialize model */
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
    int *s[2], sl[2], who;
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
  ab->dof = ab->dof0 = (ab->d == 2) ? (ab->n - 2) : (2*ab->n - 5);
  if (verbose) {
    printf("n = %3d, d = %d, dof = %3d: ", ab->n, ab->d, ab->dof);
    for (i = 0; i < ab->n; i++)
      printf("%c", ab->type[i]+'A');
    printf("\n");
  }

  nd = ab->n * ab->d;
  xnew(ab->x, nd);
  xnew(ab->x1, nd);
  xnew(ab->v, nd);
  xnew(ab->f, nd);
  xnew(ab->dx, nd);
  xnew(ab->lmx, nd);
  xnew(ab->xmin, nd);

  xnew(ab->xx[0], nd * AB_XXCNT);
  for (i = 1; i < AB_XXCNT; i++)
    ab->xx[i] = ab->xx[0] + i * nd;

#ifdef _OPENMP
  {
    int sz, ir, jr, k;
    abpairid_t *pr;
    
    ab->nthreads = omp_get_max_threads();
    xnew(ab->f_l, nd * ab->nthreads);
  
    /* partition home atoms, for thread i: [ab->homeid[i], ab->homeid[i+1]) */
    xnew(ab->homeid, ab->nthreads + 1);
    sz = (ab->n + ab->nthreads - 1) / ab->nthreads; /* chunck size */
    for (i = 0; i < ab->nthreads; i++)
      ab->homeid[ i ] = sz * i;
    ab->homeid[ ab->nthreads ] = ab->n;
 
    /* make a list of pairs */
    xnew(ab->pair, ab->n * (ab->n - 1) / 2 * sizeof(*ab->pair));
    pr = ab->pair;
    for (i = 0; i < ab->n - 2; i++) {
      for (ir = 0; ir < ab->nthreads; ir++) /* home thread */
        if (i >= ab->homeid[ ir ] && i < ab->homeid[ ir+1 ]) break;
      for (j = i + 2; j < ab->n; j++) {
        for (jr = ir; jr < ab->nthreads; jr++) /* home thread */
          if (j >= ab->homeid[ jr ] && j < ab->homeid[ jr+1 ]) break;
        pr->i = i;
        pr->j = j;
        pr->c = ab->clj[ ab->type[i] ][ ab->type[j] ];
        pr->tid =  (jr > ir && (i + j) % 2 == 0) ? jr : ir;
        pr++;
      }
    }
    ab->paircnt = pr - ab->pair;
    qsort(ab->pair, ab->paircnt, sizeof(*ab->pair), &ab_prcmp);

    /* set up pairid, [pairid[tid], pairid[tid + 1]) pairs for thread tid */
    xnew(ab->pairid, ab->nthreads + 1);
    for (k = 0, ir = 0; ir < ab->nthreads; ir++) {
      for (; k < ab->paircnt; k++)
        if (ab->pair[k].tid == ir) break;
      ab->pairid[ir] = k;
    }
    ab->pairid[ ab->nthreads ] = ab->paircnt;
  }
#endif
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
  ab_shiftcom(ab, x);
  die_if (ab_checkconn(ab, x, 0) != 0, "initpos failed, with del = %g\n", del);
  return 0;
}

/* close ab */
void ab_close(abpro_t *ab)
{
  if (!ab) return;
  free(ab->type);
  free(ab->x);
  free(ab->x1);
  free(ab->dx);
  free(ab->v);
  free(ab->f);
  free(ab->lmx);
  free(ab->xmin);
  free(ab->xx[0]);
#ifdef _OPENMP
  free(ab->f_l);
#endif
  free(ab);
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

/* shift center of x to the origin, 
 * remove center velocity and angular momentum */
void ab_rmcom(abpro_t *ab, real *x, real *v)
{
  ab_shiftcom(ab, x);
  ab_shiftcom(ab, v);
  ab_shiftang(ab, x, v); /* remove angular momentum */
}

/* write position file (which may include velocity) */
int ab_writepos(abpro_t *ab, const real *x, const real *v, const char *fname)
{
  FILE *fp;
  int i, j, d = ab->d, n = ab->n;

  if (fname == NULL) fname = "ab.pos";
  xfopen(fp, fname, "w", return -1);

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
  xfopen(fp, fname, "r", return -1);

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
    if (fgets(s, sizeof s, fp) == NULL) goto ERR;
    if (strlen(s) < 10) goto ERR;
    for (p = s, j = 0; j < d; j++, p += next) {
      if (1 != sscanf(p, fmt, x+i*d+j, &next)) {
        fprintf(stderr, "cannot read i = %d, j = %d\n", i, j);
        goto ERR;
      }
    }
    if (hasv) {
      vi = (v != NULL) ? (v + i*d) : vtmp;
      for (j = 0; j < d; j++, p += next) {
        if (1 != sscanf(p, fmt, vi+j, &next)) {
          fprintf(stderr, "cannot read i = %d, j = %d\n", i, j);
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

/* shake with additional constraints */
static int ab_shake3d(abpro_t *ab, crv3_t *x0, rv3_t *x1, rv3_t *v, real dt,
    int itmax, double tol, int verbose)
{
  int i, j, k, again, it, n = ab->n, lgcon = ab->lgcon, lgcnt = ab->lgcnt;
  real dxi[3], g, r2, r2bad, r2ref, tmp, *lgdx0;
  rv3_t *dx0 = (rv3_t *) ab->dx;
  static const real glow = .5, r2max = 4.0;
  lgconstr_t *lgc = ab->lgc;

#pragma omp threadprivate(r2max, glow)

  /* pre-compute reference difference */
  for (i = 0; i < n-1; i++)
    rv3_diff(dx0[i], x0[i+1], x0[i]);

  /* compute distance in x0 for local constraints */
  for (k = 0; lgcon && k < lgcnt; k++) {
    if (!lgc[k].on) continue;
    i = lgc[k].i;
    j = lgc[k].j;
    rv3_diff(lgc[k].dx0, x0[j], x0[i]);
  }

  for (it = 0; it < itmax; it++) {
    again = 0;
#pragma omp parallel firstprivate(n) private(i, r2, dxi, g)
   {
#ifdef _OPENMP
    int ip = omp_get_thread_num();
    int np = omp_get_num_threads();
    int n1 = n - 1;
    int sz = (n1 + np - 1)/np;
    int imin = ip * sz, imax = (ip + 1) * sz;
    if (imax > n1) imax = n1;
#else
    int imin = 0, imax = n - 1;
#endif

    for (i = imin; i < imax; i++) { /* standard constaints */
      r2 = rv3_sqr(rv3_diff(dxi, x1[i+1], x1[i]));
      if (r2 > r2max) { /* too large, impossible to correct */
        if (verbose)
          fprintf(stderr, "shake: r(%d, %d) = %g\n", i,i+1, sqrt(r2));
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
        g = (1 - r2) / (4 * g);
        rv3_sinc(x1[i],   dx0[i], -g);
        if (v) rv3_sinc(v[i], dx0[i], -g/dt);
        if (i == imax - 1) {
#pragma omp critical
         {
          rv3_sinc(x1[i+1], dx0[i],  g);
          if (v) rv3_sinc(v[i+1], dx0[i],  g/dt);
         }
        } else {
          rv3_sinc(x1[i+1], dx0[i],  g);
          if (v) rv3_sinc(v[i+1], dx0[i],  g/dt);
        }
      }
    }
   } /* end of parallel code */
    
    /* local constraints */
    if (lgcon) {
      for (k = 0; k < lgcnt; k++) {
        if (!lgc[k].on) continue;
        i = lgc[k].i;
        j = lgc[k].j;
        r2 = rv3_sqr(rv3_diff(dxi, x1[j], x1[i]));
        tmp = (r2ref = lgc[k].r2ref) * r2max;
        if (r2 > tmp) r2 = tmp;
        if (fabs(r2 - r2ref) > tol * r2ref) {
          if (!again) again = 1;
          g = rv3_dot(dxi, lgc[k].dx0);
          tmp = glow * r2ref;
          if (fabs(g) < tmp) g = g > 0 ? tmp : -tmp;
          g = (r2ref - r2) / (4 * g);
          lgdx0 = lgc[k].dx0;
          rv3_sinc(x1[i], lgdx0, -g);
          rv3_sinc(x1[j], lgdx0, +g);
          if (v) {
            rv3_sinc(v[i], lgdx0, -g/dt);
            rv3_sinc(v[j], lgdx0, +g/dt);
          }
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
  int i, j, k, again, it, n = ab->n, lgcon = ab->lgcon, lgcnt = ab->lgcnt;
  real dv[3], g, rvbad, *lgdx0;
  rv3_t *dx = (rv3_t *) ab->dx;
  lgconstr_t *lgc = ab->lgc;

  for (i = 0; i < n-1; i++)
    rv3_diff(dx[i], x0[i+1], x0[i]);

  for (k = 0; lgcon && k < lgcnt; k++) { /* local constraints */
    if (!lgc[k].on) continue;
    i = lgc[k].i;
    j = lgc[k].j;
    rv3_diff(lgc[k].dx0, x0[j], x0[i]);
  }

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

    for (k = 0; lgcon && k < lgcnt; k++) { /* local constraints */
      if (!lgc[k].on) continue;
      i = lgc[k].i;
      j = lgc[k].j;
      rv3_diff(dv, v[j], v[i]);
      lgdx0 = lgc[k].dx0;
      g = .5f * rv3_dot(lgdx0, dv);
      if (fabs(g) > tol) {
        if (!again) again = 1;
        rv3_sinc(v[i], lgdx0, +g);
        rv3_sinc(v[j], lgdx0, -g);
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
  rv3_t *dx0 = (rv3_t *) ab->xx[0], *dx1 = (rv3_t *) ab->xx[1], *x = (rv3_t *) ab->xx[2];
  real *dl = ab->xx[3], *dm = dl + n, *du = dm + n, *lam = du + n, *rhs = lam + n;
  real y;

  nl = n - 1;
  for (i = 0; i < nl; i++) {
    rv3_diff(dx1[i], x1[i], x1[i+1]);
    rv3_diff(dx0[i], x0[i], x0[i+1]);
    rhs[i] = 1 - rv3_sqr(dx1[i]);
  }

  /* dm[0..nl-1], du[0..nl-2], dl[1..nl-1] */
  dm[0] =  4*rv3_dot(dx0[0], dx1[0]);
  du[0] = -2*rv3_dot(dx0[1], dx1[0]);
  for (i = 1; i < nl; i++) {
    dl[i] = -2*rv3_dot(dx1[i], dx0[i-1]);
    dm[i] =  4*rv3_dot(dx1[i], dx0[i]);
    du[i] = -2*rv3_dot(dx1[i], dx0[i+1]); /* no dx0[nl], but doesn't matter */ 
  }

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
  rv3_t *dx = (rv3_t *) ab->xx[0], *dv = (rv3_t *) ab->xx[1];
  real *dl = ab->xx[2], *dm = dl + n, *du = dm + n, *lam = du + n, *rhs = lam + n;

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
  real ua = 0, ulj = 0;
  rv3_t *dx = (rv3_t *) ab->dx;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 0; i < n - 2; i++)
    ua += 1.f - rv3_dot(dx[i+1], dx[i]);

#pragma omp parallel for reduction(+:ulj) private(i, j)
  for (i = 0; i < n - 2; i++) {
    real dr, dr2, dr6;
    for (j = i+2; j < n; j++) {
      dr2 = rv3_dist2(r[j], r[i]);
      if (soft && dr2 < 1.f) {
        dr = (real) sqrt(dr2);
        ulj += (52 - 48*dr) - ab->clj[ab->type[i]][ab->type[j]]*(28 - 24*dr);
      } else {
        dr2 = 1/dr2;
        dr6 = dr2*dr2*dr2;
        ulj += 4*dr6*(dr6 - ab->clj[ab->type[i]][ab->type[j]]);
      }
    }
  }
  return ua * .25f  + ulj;
}

static real ab_energy3dm2(abpro_t *ab, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  real ua = 0, ud = 0, ulj = 0;
  rv3_t *dx = (rv3_t *) ab->dx;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 1; i < n-1; i++)
    ua += rv3_dot(dx[i], dx[i-1]);

  for (i = 1; i < n-2; i++)
    ud -= .5f * rv3_dot(dx[i+1], dx[i-1]);

#pragma omp parallel for reduction(+:ulj) private(i, j)
  for (i = 0; i < n-2; i++) {
    real dr2, dr6;
    for (j = i+2; j < n; j++) {
      dr2 = rv3_dist2(r[j], r[i]);
      if (soft && dr2 < 1.f) {
        ulj += ab->clj[ab->type[i]][ab->type[j]]*(ab->sla - ab->slb*(real)sqrt(dr2));
      } else {
        dr2 = 1.f/dr2;
        dr6 = dr2*dr2*dr2;
        ulj += 4.f*ab->clj[ab->type[i]][ab->type[j]]*dr6*(dr6 - 1.f);
      }
    }
  }
  return ua + ud + ulj;
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

static real ab_force3dm1(abpro_t *ab, rv3_t *f_g, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  rv3_t *dx = (rv3_t *) ab->dx;
  real U = 0.0f, ua = 0.0f, *dxm, *dxp;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 0; i < n; i++) rv3_zero(f_g[i]);
  for (i = 1; i < n - 1; i++) {
    dxp = dx[i];
    dxm = dx[i-1];
    ua += 1.f - rv3_dot(dxp, dxm);
    rv3_sinc(f_g[i-1], dxp, -.25f);
    rv3_sinc(f_g[i],   dxp,  .25f);
    rv3_sinc(f_g[i],   dxm, -.25f);
    rv3_sinc(f_g[i+1], dxm,  .25f);
  }

#pragma omp parallel firstprivate(n) private(i, j)
{ /* parallel code starts here */
  real dr2, dr6, ff, c;
  rv3_t *f, dxi;
  real ulj = 0.f;

#ifdef _OPENMP
  int ip = omp_get_thread_num();
  int np = omp_get_num_threads();
  int imin = ab->homeid[ip], imax = ab->homeid[ip + 1];
  int ipr, ipr0 = ab->pairid[ip], ipr1 = ab->pairid[ip + 1];
  rv3_t *f_l = (rv3_t *) ab->f_l;
  f = f_l + n * ip; /* point to the proper local force */
  /* clear the local force */
  for (i = 0; i < n; i++) rv3_zero(f[i]);
#else
  f = f_g;
#endif
  
#ifdef _OPENMP
  for (ipr = ipr0; ipr < ipr1; ipr++) { 
    i = ab->pair[ipr].i;
    j = ab->pair[ipr].j;
    c = ab->pair[ipr].c;
#else
  for (i = 0; i < n - 2; i++) 
  for (j = i + 2; j < n; j++) {
    c = ab->clj[ ab->type[i] ][ ab->type[j] ];
#endif
  
    dr2 = rv3_sqr( rv3_diff(dxi, r[i], r[j]) );

    if (soft && dr2 < 1.) {
      dr2 = (real) sqrt(dr2);
      ulj += (52.f - 28.f*c) - 24.f*dr2*(2.f-c);
      ff = 24.f*(2.f - c)/dr2;
    } else {
      dr2 = 1.f/dr2;
      dr6 = dr2*dr2*dr2;
      ulj += 4.f*dr6*(dr6 - c);
      ff = 24.f*dr2*dr6*(dr6*2.f - c);
    }
    rv3_sinc(f[i], dxi, +ff);
    rv3_sinc(f[j], dxi, -ff);
  }
#ifdef _OPENMP
#pragma omp barrier
  for (i = imin; i < imax; i++) /* collect global force, f_l flushed by barrier */
    for (j = 0; j < np; j++) { 
      rv3_inc(f_g[i], f_l[n*j + i]);
    }
#endif
#pragma omp atomic 
  U += ulj;
} /* parallel code stops */
  return ua * 0.25f + U;
}

static real ab_force3dm2(abpro_t *ab, rv3_t *f_g, crv3_t *r, int soft)
{
  int i, j, n = ab->n;
  rv3_t *dx = (rv3_t *) ab->dx;
  real U = 0.f, ua = 0.f, ud = 0.f, *dxm, *dxp;

  for (i = 0; i < n - 1; i++)
    rv3_diff(dx[i], r[i+1], r[i]);

  for (i = 0; i < n; i++) rv3_zero(f_g[i]);

  for (i = 1; i < n-1; i++) {
    dxp = dx[i];
    dxm = dx[i-1];
    rv3_inc(f_g[i-1], dxp);
    rv3_dec(f_g[i],   dxp);
    rv3_inc(f_g[i],   dxm);
    rv3_dec(f_g[i+1], dxm);
    ua += rv3_dot(dxp, dxm);
  }
  
  for (i = 1; i < n-2; i++) {
    dxp = dx[i+1];
    dxm = dx[i-1];
    rv3_sinc(f_g[i-1], dxp, -.5f);
    rv3_sinc(f_g[i],   dxp,  .5f);
    rv3_sinc(f_g[i+1], dxm, -.5f);
    rv3_sinc(f_g[i+2], dxm,  .5f);
    ud -= .5f*rv3_dot(dxp, dxm);
  }

#pragma omp parallel firstprivate(n) private(i, j) 
{ /* parallel code starts here */
  real dr2, dr6, ff, c;
  rv3_t *f, dxi;
  real ulj = 0.f;

#ifdef _OPENMP
  int ip = omp_get_thread_num();
  int np = omp_get_num_threads();
  int imin = ab->homeid[ip], imax = ab->homeid[ip + 1];
  int ipr, ipr0 = ab->pairid[ip], ipr1 = ab->pairid[ip + 1];
  rv3_t *f_l = (rv3_t *) ab->f_l;
  f = f_l + n * ip; /* point to the proper local force */
  /* clear the local force */
  for (i = 0; i < n; i++) rv3_zero(f[i]);
#else
  f = f_g;
#endif

#ifdef _OPENMP
  for (ipr = ipr0; ipr < ipr1; ipr++) { 
    i = ab->pair[ipr].i;
    j = ab->pair[ipr].j;
    c = ab->pair[ipr].c;
#else
  for (i = 0; i < n - 2; i++) 
  for (j = i + 2; j < n; j++) {
    c = ab->clj[ ab->type[i] ][ ab->type[j] ];
#endif

    dr2 = rv3_sqr( rv3_diff(dxi, r[i], r[j]) );

    if (soft && dr2 < 1.f) {
      dr2 = (real) sqrt(dr2);
      ulj += c*(ab->sla - ab->slb*dr2);
      ff = (ab->slb*c)/dr2;
    } else {
      dr2 = 1.f/dr2;
      dr6 = dr2*dr2*dr2;
      ulj += 4.f*c*dr6*(dr6 - 1.f);
      ff = 48.f*c*dr2*dr6*(dr6 - .5f);
    }
    rv3_sinc(f[i], dxi, +ff);
    rv3_sinc(f[j], dxi, -ff);
  }
#ifdef _OPENMP
#pragma omp barrier
  for (i = imin; i < imax; i++) /* collect global force, barrier flushes f_l */
    for (j = 0; j < np; j++) 
      rv3_inc(f_g[i], f_l[n*j + i]);
#endif
#pragma omp atomic 
  U += ulj;
}
  return ua + ud + U;
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
  real **x = ab->xx, **f = ab->xx + 2, *v = ab->xx[4];
  const real DELMAX = 0.20f;

  if (itmax <= 0) itmax = 10000;
  if (tol <= 0.) tol = 1e-12;
  /* to make a working copy */
  memcpy(x[id = 0], r, n*d*sizeof(real));
  up = ab_force(ab, f[id], x[id], 0);
  memset(v, 0, n*d*sizeof(real));

  for (t = 1; t <= itmax; t++) {
    for (i = 0; i < n; i++)
      for (j = 0; j < d; j++) {
        del = v[i*d+j] = v[i*d+j] * mem + f[id][i*d+j] * step;
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

static int ab_vv3d(abpro_t *ab, real fscal, real dt, int soft, int milc)
{
  int i, verbose = 1, n = ab->n;
  real dth = .5f*dt*fscal;
  rv3_t *v = (rv3_t *)ab->v, *x = (rv3_t *)ab->x, *x1 = (rv3_t *)ab->x1, *f = (rv3_t *)ab->f;

#pragma omp parallel for schedule(static)
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
  
#pragma omp parallel for schedule(static)
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

/* local geometric constraints (LGC) */

/* initialize local constraints,
 * level 1: aba, 2: abba, 3: both */
INLINE void ab_initconstr(abpro_t *ab, int level)
{
  int i, di, n = ab->n;
  real dr0, rref[2][2] = {{1.115f, 1.16f}, {1.09f, 1.2f}}; /* empirical constraint distances */

  if (level < 0) {
    level = 3; /* turn on all optimizations */
    if (ab->model == 2) level &= ~2;
  }
  if (ab->d == 2 || level == 0) return;
  ab->lgcnt = 0;
  xnew(ab->lgc, 1);
  for (di = 2; di <= 3; di++) { /* loop over ABA and ABBA */
    if ( !(level & (di - 1)) ) continue;
    dr0 = rref[ab->model - 1][di - 2];

    for (i = 0; i < n - di; i++) {
      if (ab->type[i] != 0 || ab->type[i + di] != 0) continue;
      xrenew(ab->lgc, ab->lgcnt + 1);
      ab->lgc[ab->lgcnt].i = i;
      ab->lgc[ab->lgcnt].j = i + di;
      ab->lgc[ab->lgcnt].on = 0;
      ab->lgc[ab->lgcnt].r2ref = dr0 * dr0;
      ab->lgcnt++;
    }
  }
  ab->lgcon = 1; /* turn on constraints by default */
  ab->lgact = 0; /* no constraint is active yet */
}

/* update constraints if atoms are close enough*/
INLINE void ab_updconstr(abpro_t *ab)
{
  int i, j, k, lgcnt = ab->lgcnt;
  real dr2;

  if (!ab->lgcon) return;
  for (k = 0; k < lgcnt; k++) { /* try to turn on constraints */
    if (ab->lgc[k].on) continue;
    i = ab->lgc[k].i;
    j = ab->lgc[k].j;
    dr2 = rv3_dist2(ab->x + 3*i, ab->x + 3*j);
    if (fabs(dr2 - ab->lgc[k].r2ref) < 0.3) {
      ab->lgc[k].on = 1;
      ab->lgact++;
    }
  }
  ab->dof = ab->dof0 - ab->lgact;
}
#endif

