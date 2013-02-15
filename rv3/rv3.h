#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef RESTRICT
#define RESTRICT __restrict 
#endif
#include "def.h"
#include "rng.c"
#ifndef RV3_H__
#define RV3_H__

#ifndef DV3_T
#define DV3_T dv3_t
  typedef double dv3_t[3];
  typedef const double cdv3_t[3];
  typedef double dm3_t[3][3];
#endif

#ifndef RV3_T
#define RV3_T rv3_t
  typedef real rv3_t[3];
  typedef const real crv3_t[3];
  typedef real rm3_t[3][3];
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define rv3_fprint_(fp, r, nm, fmt, nl) { int i_; \
  if (nm) fprintf(fp, "%s: ", nm); \
  for (i_ = 0; i_ < 3; i_++) fprintf(fp, fmt, r[i_]); \
  fprintf(fp, "%c", (nl ? '\n' : ';')); }

#define rm3_fprint_(fp, m, nm, fmt, nl) { int i_, j_; \
  if (nm) fprintf(fp, "%s:%c", nm, (nl ? '\n' : ' ')); \
  for (i_ = 0; i_ < 3; i_++) { \
    for (j_ = 0; j_ < 3; j_++) \
      fprintf(fp, fmt, m[i_][j_]); \
    fprintf(fp, "%s", (nl ? "\n" : "; ")); } }

#define dv3_print(r, nm, fmt, nl) dv3_fprint(stdout, r, nm, fmt, nl)
INLINE void dv3_fprint(FILE *fp, double *r, const char *nm, const char *fmt, int nl)
  rv3_fprint_(fp, r, nm, fmt, nl)

#define dm3_print(m, nm, fmt, nl) dm3_fprint(stdout, m, nm, fmt, nl)
INLINE void dm3_fprint(FILE *fp, double (*m)[3], const char *nm, const char *fmt, int nl)
  rm3_fprint_(fp, m, nm, fmt, nl)

#define rv3_print(r, nm, fmt, nl) rv3_fprint(stdout, r, nm, fmt, nl)
INLINE void rv3_fprint(FILE *fp, real *r, const char *nm, const char *fmt, int nl)
  rv3_fprint_(fp, r, nm, fmt, nl)

#define rm3_print(m, nm, fmt, nl) rm3_fprint(stdout, m, nm, fmt, nl)
INLINE void rm3_fprint(FILE *fp, real (*m)[3], const char *nm, const char *fmt, int nl)
  rm3_fprint_(fp, m, nm, fmt, nl)

/* due to possible pointer overlap, be careful when to use 'const' */

INLINE double *dv3_make(double *x, double a, double b, double c)
  { x[0] = a; x[1] = b; x[2] = c; return x; }
INLINE double *dv3_fromrv3(double *x, const real *rx)
  { return dv3_make(x, rx[0], rx[1], rx[2]); }
INLINE double *dv3_zero(double *x) { return dv3_make(x, 0, 0, 0); }
INLINE double *dv3_copy(double *x, const double *src)
  { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; return x; }
#define dv3_ncopy(x, src, n) memcpy(x, src, 3*n*sizeof(double))
INLINE void dv3_swap(double *x, double *y)
  { double z[3]; dv3_copy(z, x); dv3_copy(x, y); dv3_copy(y, z); }

INLINE double dv3_sqr (const double *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
INLINE double dv3_norm(const double *x) { return (double)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

INLINE double dv3_dot(const double *x, const double *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

INLINE double *dv3_cross(double *RESTRICT z, const double *x, const double *y)
{
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
  return z;
}

INLINE double *dv3_neg(double *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  x[2] = -x[2];
  return x;
}

INLINE double *dv3_neg2(double *nx, const double *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  nx[2] = -x[2];
  return nx;
}

INLINE double *dv3_inc(double * RESTRICT x, const double *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  x[2] += dx[2];
  return x;
}

INLINE double *dv3_dec(double *x, const double *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  x[2] -= dx[2];
  return x;
}

INLINE double *dv3_sinc(double * RESTRICT x, const double *dx, double s)
{
  x[0] += dx[0]*s;
  x[1] += dx[1]*s;
  x[2] += dx[2]*s;
  return x;
}

INLINE double *dv3_smul(double *x, double s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}

INLINE double *dv3_smul2(double * RESTRICT y, const double *x, double s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  y[2] = x[2]*s;
  return y;
}

INLINE double *dv3_normalize(double *x)
{
  double r = dv3_norm(x);
  if (r > 0.0) dv3_smul(x, 1.0/r);
  return x;
}

INLINE double *dv3_makenorm(double *v, double x, double y, double z)
  { return dv3_normalize( dv3_make(v, x, y, z) ); }

INLINE double *dv3_diff(double * RESTRICT diff, const double *a, const double *b)
{
  diff[0] = a[0] - b[0];
  diff[1] = a[1] - b[1];
  diff[2] = a[2] - b[2];
  return diff;
}

INLINE double dv3_dist2(const double *a, const double *b) 
{
  double d[3]; 
  return dv3_sqr(dv3_diff(d, a, b));
}

INLINE double dv3_dist(const double *a, const double *b) 
{
  return (double) sqrt(dv3_dist2(a, b));
}

INLINE double *dv3_add(double * RESTRICT sum, const double *a, const double *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  sum[2] = a[2]+b[2];
  return sum;
}

INLINE double *dv3_nadd(double *sum, const double *a, const double *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  sum[2] = -a[2]-b[2];
  return sum;
}

INLINE double *dv3_lincomb2(double * RESTRICT sum, const double *a, const double *b, double s1, double s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  sum[2] = a[2]*s1+b[2]*s2;
  return sum;
}

INLINE real *rv3_make(real *x, real a, real b, real c)
  { x[0] = a; x[1] = b; x[2] = c; return x; }
INLINE real *rv3_fromdv3(real *x, const double *dx)
  { return rv3_make(x, (real) dx[0], (real) dx[1], (real) dx[2]); }
INLINE real *rv3_zero(real *x) { return rv3_make(x, 0, 0, 0); }
INLINE real *rv3_copy(real *x, const real *src)
  { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; return x; }
/* use macro to avoid const qualifier of src */
#define rv3_ncopy(x, src, n) memcpy(x, src, 3*n*sizeof(real))
INLINE void rv3_swap(real *x, real *y)
  { real z[3]; rv3_copy(z, x); rv3_copy(x, y); rv3_copy(y, z); }

INLINE real rv3_sqr (const real *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
INLINE real rv3_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

/* if x == y, try to use sqr */
INLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

INLINE real *rv3_cross(real *RESTRICT z, const real *x, const real *y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}

INLINE real *rv3_neg(real *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  x[2] = -x[2];
  return x;
}

INLINE real *rv3_neg2(real *nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  nx[2] = -x[2];
  return nx;
}

INLINE real *rv3_inc(real * RESTRICT x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  x[2] += dx[2];
  return x;
}

INLINE real *rv3_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  x[2] -= dx[2];
  return x;
}

INLINE real *rv3_sinc(real * RESTRICT x, const real *dx, real s)
{
  x[0] += dx[0]*s;
  x[1] += dx[1]*s;
  x[2] += dx[2]*s;
  return x;
}

INLINE real *rv3_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}

/* if y == x, just use smul */
INLINE real *rv3_smul2(real * RESTRICT y, const real *x, real s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  y[2] = x[2]*s;
  return y;
}

INLINE real *rv3_normalize(real *x)
{
  real r = rv3_norm(x);
  if (r > 0.f) rv3_smul(x, 1.f/r);
  return x;
}

INLINE real *rv3_makenorm(real *v, real x, real y, real z)
  { return rv3_normalize( rv3_make(v, x, y, z) ); }

/* for in-place difference use rv3_dec */
INLINE real *rv3_diff(real * RESTRICT diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  diff[2] = a[2]-b[2];
  return diff;
}

/* distance^2 between a and b */
INLINE real rv3_dist2(const real *a, const real *b) 
{
  real d[3]; 
  return rv3_sqr(rv3_diff(d, a, b));
}

/* distance between a and b */
INLINE real rv3_dist(const real *a, const real *b) 
{
  return (real) sqrt(rv3_dist2(a, b));
}

/* sum = a+b, for in-place addition use rv3_inc */
INLINE real *rv3_add(real * RESTRICT sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  sum[2] = a[2]+b[2];
  return sum;
}

/* sum = -a-b */
INLINE real *rv3_nadd(real *sum, const real *a, const real *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  sum[2] = -a[2]-b[2];
  return sum;
}

INLINE real *rv3_lincomb2(real * RESTRICT sum, const real *a, const real *b, real s1, real s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  sum[2] = a[2]*s1+b[2]*s2;
  return sum;
}

/* angle and gradients of cos(x1-x2-x3) */
INLINE real rv3_cosang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real a[3], b[3], ra, rb, dot;

  ra = rv3_norm(rv3_diff(a, x1, x2));
  rv3_smul(a, 1.f/ra);
  rb = rv3_norm(rv3_diff(b, x3, x2));
  rv3_smul(b, 1.f/rb);
  dot = rv3_dot(a, b);
  if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
  if (g1) {
    rv3_lincomb2(g1, b, a, 1.f/ra, -dot/ra);
    rv3_lincomb2(g3, a, b, 1.f/rb, -dot/rb);
    rv3_nadd(g2, g1, g3);
  }
  return dot;
}

/* angle and gradients of x1-x2-x3 */
INLINE real rv3_ang(const real *x1, const real *x2, const real *x3,
    real *g1, real *g2, real *g3)
{
  real dot, sn;
  
  dot = rv3_cosang(x1, x2, x3, g1, g2, g3);
  sn = (real) sqrt(1 - dot*dot);
  if (sn < 1e-7) sn = 1; else sn = -1.f/sn;
  if (g1) {
    rv3_smul(g1, sn);
    rv3_smul(g2, sn);
    rv3_smul(g3, sn);
  }
  return (real) acos(dot);
}

/* vertical distance from x to line a-b */
INLINE real rv3_vdist(const real *x, const real *a, const real *b)
{
  real nm[3], d[3], dot;

  rv3_diff(d, x, a);
  rv3_normalize(rv3_diff(nm, a, b));
  dot = rv3_dot(d, nm);
  return rv3_norm(rv3_sinc(d, nm, -dot));
}

/* signed distance from x to the plane extended by a, b, c */
INLINE real rv3_vpdist(const real *x, const real *a, const real *b, const real *c)
{
  real u[3], v[3], m[3];

  rv3_diff(u, b, a);
  rv3_diff(v, c, b);
  rv3_normalize(rv3_cross(m, u, v));
  rv3_diff(u, x, a);
  return rv3_dot(u, m);
}

/* light weight dihedral */
INLINE real rv3_dih(const real xi[], const real xj[], const real xk[], const real xl[],
    real gi[], real gj[], real gk[], real gl[])
{
  real tol, phi, cosphi = 1.f;
  real nxkj, nxkj2, m2, n2;
  real xij[3], xkj[3], xkl[3], uvec[3], vvec[3], svec[3];
  real m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */

  rv3_diff(xij, xi, xj);
  rv3_diff(xkj, xk, xj);
  rv3_diff(xkl, xk, xl);
  nxkj2 = rv3_sqr(xkj);
  nxkj = (real) sqrt(nxkj2);
  tol = (sizeof(real) == sizeof(float)) ? nxkj2 * 6e-8f : nxkj2 * 1e-16f;

  rv3_cross(m, xij, xkj);
  m2 = rv3_sqr(m);
  rv3_cross(n, xkj, xkl);
  n2 = rv3_sqr(n);
  if (m2 > tol && n2 > tol) {
    cosphi = rv3_dot(m, n);
    cosphi /= (real) sqrt(m2 * n2);
    if (cosphi >= 1.f) cosphi = 1.f;
    else if (cosphi < -1.f) cosphi = -1.f;
  }
  phi = (real) acos(cosphi);
  if (rv3_dot(n, xij) < 0.0) phi = -phi;

  /* optionally calculate the gradient */
  if (gi != NULL) {
    if (m2 > tol && n2 > tol) {
      rv3_smul2(gi, m, nxkj/m2);
      rv3_smul2(gl, n, -nxkj/n2);
      rv3_smul2(uvec, gi, rv3_dot(xij, xkj)/nxkj2);
      rv3_smul2(vvec, gl, rv3_dot(xkl, xkj)/nxkj2);
      rv3_diff(svec, uvec, vvec);
      rv3_diff(gj, svec, gi);
      rv3_nadd(gk, svec, gl);
    } else { /* clear the gradients */
      rv3_zero(gi);
      rv3_zero(gj);
      rv3_zero(gk);
      rv3_zero(gl);
    }
  }
  return phi;
}


/* structure for dihedral calculation */
typedef struct {
  int  szreal; /* sizeof real */
  int  pad0;   /* padding */
  real phi; /* cis is zero, clockwise positive */
  real cos; /* cos(m, n) */
  real sgn; /* (0, pi) is 1.0, otherwise -1.0 */

  real g2;
  real g[4][3]; /* gradient for each particle */

  real div; /* the divengence */
  real d4ij, d4ik, d4jj, d4jk, d4jl, d4kk, d4kl;

  unsigned int flags; /* a copy of flags used */
  int t1, t2, t3; /* gromacs shift indices */
  const void *pbcdata; /* periodic boundary condition descriptor */
  int (*pbcdiff)(real *xij, const real *xi, const real *xj, const void *); 
    /* function to handle pbc, use GROMACS convention: the last is the difference */
} dihcalc_t;

#define DIH_GRAD  0x0001
#define DIH_DIV   0x0002
/*#define DIH_CONJ  0x0004 */
/*#define DIH_PROJ  0x0008 */
#define DIH_I     0x0010
#define DIH_J     0x0020
#define DIH_K     0x0040
#define DIH_L     0x0080
#define DIH_FOUR  (DIH_I|DIH_J|DIH_K|DIH_L)
/* the four atoms involved */
#define DIH_ALL   (DIH_FOUR|DIH_GRAD|DIH_DIV)
/* only I and L, so no divergence */
#define DIH_ENDS  (DIH_GRAD|DIH_I|DIH_L)
/* polymer convention, 0 == trans */
#define DIH_POLYMER 0x1000

/* compute the dihedral angle, gradient g and divegence
 * of the field v conjugate to gradient (v.g = 1)
 *
 * if dih is NULL and flags is 0, only the dihedral angle is computed
 * optionally, the gradient and divergent are computed with flags
 * DIH_GRAD and DIH_DIV respectively (the latter requires the former)
 * routine for treating periodic boundary condition can specified by
 * assigning a function pointer to dih->pbcdiff with additional information
 * to dih->pbcdata entry, otherwise, dih->pbcdiff *must* be set to NULL
 * the procedure of computing similar to that in GROMACS
 *
 * the conjugate field v = g / [g.g], such that v.g = 1, where g = grad(phi)
 * it is not explicitly computed (since it lies along the gradient)
 * however, the denominator is saved to dih->g2
 * the calculation of the divergence of v is simplified by the fact that
 * the divergence of the gradient is always zero, i.e., div.g = 0, thus
 *     div.v = -2 [ g.(gg).g ] /[g.g]^2,
 * where gg = grad grad(phi) involves: d4ij, d4ik, ..., d4kl.
 *
 * both g and div can be computed for a subset of the four involving atoms
 * by passing `flags' a combination of DIH_I, DIH_J, DIH_K and DIH_L
 * however, *all* moments d4ij, d4ik, ... , d4kl are always calculated
 * though only the involved ones are used to produce the divergence. */
#define rv3_calcdihv(dih, x, idx, flags) \
  rv3_calcdih(dih, x[*(idx)], x[*(idx+1)], x[*(idx+2)], x[*(idx+3)], flags)

INLINE real rv3_calcdih(dihcalc_t *dih,
    const real *xi, const real *xj, const real *xk, const real *xl,
    unsigned int flags)
{
  real dot, scl, tol, vol, phi, sgn, cosphi;
  real nxkj, nxkj2, m2, n2;
  real xij[3], xkj[3], xkl[3];
  real m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */

  if (dih != NULL && sizeof(real) != dih->szreal) {
    fprintf(stderr, "size real should be %d instead of %d\n",
        (int) sizeof(real), (int) dih->szreal);
    exit(1);
  }
  if (dih != NULL && dih->pbcdiff != NULL) { /* handle pbc */
    dih->t1 = (*dih->pbcdiff)(xij, xi, xj, dih->pbcdata);
    dih->t2 = (*dih->pbcdiff)(xkj, xk, xj, dih->pbcdata);
    dih->t3 = (*dih->pbcdiff)(xkl, xk, xl, dih->pbcdata);
  } else {
    rv3_diff(xij, xi, xj);
    rv3_diff(xkj, xk, xj);
    rv3_diff(xkl, xk, xl);
  }
  nxkj2 = rv3_sqr(xkj);
  nxkj = (real) sqrt(nxkj2);
  if (sizeof(real) <= 4)
    tol = nxkj2 * 6e-8f;
  else
    tol = nxkj2 * 1e-16f;

  rv3_cross(m, xij, xkj);
  m2 = rv3_sqr(m);
  rv3_cross(n, xkj, xkl);
  n2 = rv3_sqr(n);
  if (m2 > tol && n2 > tol) {
    scl = (real) sqrt(m2*n2);
    dot = rv3_dot(m, n);
    cosphi = dot/scl;
    if (cosphi >= (real) 1.) cosphi = 1.f; 
    else if (cosphi < (real)(-1.)) cosphi = -1.f;
  } else {
    cosphi = 1.f;
  }
  phi = (real) acos(cosphi);
  vol = rv3_dot(n, xij);
  sgn = (vol > 0.0f) ? 1.0f : -1.0f;
  phi *= sgn;
  if (flags & DIH_POLYMER) { /* switch to polymer convention, 0 == trans */
    if (phi > 0) phi -= M_PI;
    else phi += M_PI;
  }
  if (dih != NULL) {
    dih->phi = phi;
    dih->sgn = sgn;
    dih->cos = cosphi;
    dih->flags = flags;
  }

  /* optionally calculate the gradient */
  if (dih != NULL && (flags & (DIH_GRAD|DIH_DIV))) { /* divergence implies gradient */
    /* clear divergence */
    dih->div = dih->d4ij = dih->d4ik = dih->d4jj = dih->d4jk = dih->d4jl = dih->d4kk = dih->d4kl = 0.0f;

    /* calculate the gradient of the dihedral */
    if (m2 > tol && n2 > tol) {
      real uvec[3], vvec[3], svec[3], g2all, invg2;
      unsigned doi, doj, dok, dol;

      doi = (flags & DIH_I);
      doj = (flags & DIH_J);
      dok = (flags & DIH_K);
      dol = (flags & DIH_L);

      rv3_smul2(dih->g[0], m,  nxkj/m2);
      rv3_smul2(dih->g[3], n, -nxkj/n2);

      rv3_smul2(uvec, dih->g[0], rv3_dot(xij, xkj)/nxkj2);
      rv3_smul2(vvec, dih->g[3], rv3_dot(xkl, xkj)/nxkj2);
      rv3_diff(svec, uvec, vvec);

      rv3_diff(dih->g[1], svec, dih->g[0]);
      rv3_nadd(dih->g[2], svec, dih->g[3]);

      g2all = 0.0f;
      if (doi) g2all += rv3_sqr(dih->g[0]);
      if (doj) g2all += rv3_sqr(dih->g[1]);
      if (dok) g2all += rv3_sqr(dih->g[2]);
      if (dol) g2all += rv3_sqr(dih->g[3]);
      dih->g2 = g2all;
      invg2 = 1.0f/g2all;

      if (flags & DIH_DIV) {
        real xkjv[3], nvv[3], mvv[3];
        real gjxij, gjmvv, gjxkl, gjnvv;
        real gkmvv, gknvv, gkxkl, gkxij;
        real kivkj, klvkj, ljvkj, ijvkj;
        real kikl, ijlj;
        real tmp1, tmp2;
        real sinmn;

        rv3_smul2(mvv, m, 1.0f/m2);
        rv3_smul2(nvv, n, 1.0f/n2);
        rv3_smul2(xkjv, xkj, 1.0f/nxkj);

        sinmn = vol*nxkj/(m2*n2);

        ijvkj = rv3_dot(xij, xkjv);
        kivkj = nxkj-ijvkj;
        klvkj = rv3_dot(xkl, xkjv);
        ljvkj = nxkj-klvkj;

        ijlj = ijvkj*ljvkj;
        kikl = kivkj*klvkj;

        gjxij = rv3_dot(dih->g[1], xij);
        gjxkl = rv3_dot(dih->g[1], xkl);
        gjmvv = rv3_dot(dih->g[1], mvv);
        gjnvv = rv3_dot(dih->g[1], nvv);
        gkxij = rv3_dot(dih->g[2], xij);
        gkxkl = rv3_dot(dih->g[2], xkl);
        gkmvv = rv3_dot(dih->g[2], mvv);
        gknvv = rv3_dot(dih->g[2], nvv);

        tmp1 = nxkj2*sinmn;
        tmp2 = tmp1/m2;
        dih->d4ij = kikl*tmp2;
        dih->d4ik = ijlj*tmp2;
        tmp2 = tmp1/n2;
        dih->d4jl = kikl*tmp2;
        dih->d4kl = ijlj*tmp2;

        dih->d4jj = -(gjxij*gjmvv+gjxkl*gjnvv)/nxkj
                +2.0f*(kivkj*gjmvv-klvkj*gjnvv)*(-kikl*sinmn);

        dih->d4jk = (gjxij*gkmvv+gjxkl*gknvv)/nxkj
              +(-(gjmvv*ljvkj+gkmvv*klvkj)*(ijvkj*kivkj)
                +(gjnvv*ijvkj+gknvv*kivkj)*(ljvkj*klvkj) )*sinmn;

        dih->d4kk = -(gkxkl*gknvv+gkxij*gkmvv)/nxkj
                +2.0f*(ljvkj*gknvv-ijvkj*gkmvv)*(ijlj*sinmn);

        /* summarize */
        if ((flags & DIH_FOUR) == DIH_FOUR) {
          tmp1 = dih->d4jj + dih->d4kk;
          tmp2 = dih->d4ij + dih->d4ik+dih->d4jk+dih->d4jl+dih->d4kl;
        } else {
          tmp1 = tmp2 = 0.0f;
          if (doj) { tmp1 += dih->d4jj; }
          if (dok) { tmp1 += dih->d4kk; }
          if (doi && doj) tmp2 += dih->d4ij;
          if (doi && dok) tmp2 += dih->d4ik;
          if (doj && dok) tmp2 += dih->d4jk;
          if (doj && dol) tmp2 += dih->d4jl;
          if (dok && dol) tmp2 += dih->d4kl;
        }
        dih->div = -2.0f*(tmp1+2.0f*tmp2)*(invg2*invg2);
      } /* do divengence */

    } else { /* clear the gradients */
      int j;
      for (j = 0; j < 4; j++)
        rv3_zero(dih->g[j]);
    }
  }
  return phi;
}

/* compute the trihedral angle
 * http://planetmath.org/encyclopedia/TrihedralAngle.html
 * tan(omega/2) = vol/den,
 * where
 * den = r1 (r2.r3) + r2 (r1.r3) + r3 (r1.r2) + r1 r2 r3 
 */
INLINE real rv3_solidang(const real v1[], const real v2[], const real v3[],
  real g1[], real g2[], real g3[])
{
  real vc[3];
  real r1, r2, r3;
  real ang, vol, den, v2d2, scnum, scden;
  const real eps = 1e-10;

  r1 = rv3_norm(v1);
  r2 = rv3_norm(v2);
  r3 = rv3_norm(v3);
  if (g1 != NULL) rv3_zero(g1);
  if (g2 != NULL) rv3_zero(g2);
  if (g3 != NULL) rv3_zero(g3);

  /* at least two points coincide */
  if (r1 < eps || r2 < eps || r3 < eps) return 0;
    
  /* the numerator */
  vol = rv3_dot(v3, rv3_cross(vc, v1, v2));
  
  /* the denominator */
  den = r1*rv3_dot(v3, v2) + r2*rv3_dot(v3, v1)
      + r3*(r1*r2 + rv3_dot(v1, v2));

  v2d2 = vol*vol + den*den;
  
  /* this happens if two vector are opposite
     the solid angle could be evolved from  +/-pi or 0
     but unfortunately, we don't know which one */
  if (v2d2 < eps) return 0;
  
  if (g1 != NULL && g2 != NULL && g3 != NULL) { /* compute the gradients */
    scnum =  2.f*den/v2d2;
    scden = -2.f*vol/v2d2;
    
    /* cross products */
    rv3_smul(rv3_cross(g3, v1, v2), scnum);
    rv3_smul(rv3_cross(g1, v2, v3), scnum);
    rv3_smul(rv3_cross(g2, v3, v1), scnum);
    
    /* compute the contributions to the denominator */
    rv3_lincomb2(vc, v2, v3, r3, r2);
    rv3_sinc(vc, v1, (rv3_dot(v2, v3) + r2*r3)/r1);
    rv3_sinc(g1, vc, scden);
  
    rv3_lincomb2(vc, v1, v3, r3, r1);
    rv3_sinc(vc, v2, (rv3_dot(v1, v3) + r1*r3)/r2);
    rv3_sinc(g2, vc, scden);
  
    rv3_lincomb2(vc, v1, v2, r2, r1);
    rv3_sinc(vc, v3, (rv3_dot(v1, v2) + r1*r2)/r3);
    rv3_sinc(g3, vc, scden);
  }
  
  /* calculate tan(omega/2) */
  ang = atan2(vol, den); /* 0 to pi */
  return 2*ang;
}



/* compute the Gauss integral for the two line segments, with gradients 
 *    int_i \int_j (dri X drj).rij/ rij^3, 
 * over two line segments rip - ri, rjp - rj, and rij = ri - rj
 * (-2 pi, 2 pi)
 * Note the sign is opposite to that of the dihedral */
INLINE real rv3_solidang2g(const real ri[], const real rip[], const real rj[],
  const real rjp[], real gi[], real gip[], real gj[], real gjp[])
{
  rv3_t v0, v1, v2, v3, g0, g1, g2, g3, g4, g5;
  real ang1, ang2;

  rv3_diff(v0, ri, rj);
  rv3_diff(v1, ri, rjp);
  rv3_diff(v2, rip, rj);
  rv3_diff(v3, rip, rjp);
  
  ang1 = rv3_solidang(v0, v1, v2, g0, g1, g2);
  ang2 = rv3_solidang(v2, v1, v3, g3, g4, g5);
  
  rv3_inc(rv3_inc(rv3_copy(gi,  g0), g1), g4);
  rv3_inc(rv3_inc(rv3_copy(gip, g2), g3), g5);
  rv3_neg(rv3_inc(rv3_inc(rv3_copy(gj,  g0), g2), g3));
  rv3_neg(rv3_inc(rv3_inc(rv3_copy(gjp, g1), g4), g5));
  
  return ang1 + ang2;
}


/* compute the double integral, old code
 *    \int_i \int_j (dri X drj).rij/ rij^3, 
 * over two line segments rip - ri, rjp - rj, and rij = ri - rj
 * (-2 pi, 2 pi)
 * Note the sign is opposite to that of the dihedral */
INLINE real rv3_solidang2(const real *ri, const real *rip, 
    const real *rj, const real *rjp)
{
  real v0[3], v1[3], v2[3], v3[3], vc[3];
  double r0, r1, r2, r3;
  double ang, vol, dn1, dn2, dn, tmp;

  r0 = rv3_norm(rv3_diff(v0, ri, rj));
  r1 = rv3_norm(rv3_diff(v1, ri, rjp));
  r2 = rv3_norm(rv3_diff(v2, rip, rj));
  r3 = rv3_norm(rv3_diff(v3, rip, rjp));

  /* avoid coplanar vectors */
  vol = rv3_dot(v0, rv3_cross(vc, v1, v2));
  if(fabs(vol) < 1e-28) return 0;

  /* calculate the denominator */
  tmp = r1*r2 + rv3_dot(v1, v2);
  /* http://planetmath.org/encyclopedia/TrihedralAngle.html 
   * tan(omega/2) = vol/den,
   * where
   * den = r1 (r2.r3) + r2 (r1.r3) + r3 (r1.r2) + r1 r2 r3
   * */
  dn1 = r1*rv3_dot(v0, v2) + r2*rv3_dot(v0, v1) + r0*tmp;
  dn2 = r1*rv3_dot(v3, v2) + r2*rv3_dot(v3, v1) + r3*tmp;

  /* calculate tan(omega1/2 + omega2/2) */
  dn = (dn1 + dn2)/(dn1*dn2 - vol*vol);
  ang = atan(fabs(vol) * dn) + (dn < 0 ? M_PI : 0); /* 0 to pi */

  return (real) (vol > 0 ? 2*ang : -2*ang);
}



INLINE dv3_t *dm3_fromrm3(double a[3][3], real ra[3][3])
{
  a[0][0] = ra[0][0];
  a[0][1] = ra[0][1];
  a[0][2] = ra[0][2];
  a[1][0] = ra[1][0];
  a[1][1] = ra[1][1];
  a[1][2] = ra[1][2];
  a[2][0] = ra[2][0];
  a[2][1] = ra[2][1];
  a[2][2] = ra[2][2];
  return a;
}


/* a = b */
INLINE dv3_t *dm3_copy(double a[3][3], double b[3][3])
{
  dv3_copy(a[0], b[0]);
  dv3_copy(a[1], b[1]);
  dv3_copy(a[2], b[2]);
  return a;
}

/* transpose */
INLINE dv3_t *dm3_trans(double a[3][3]) 
{
  double x;
  x = a[0][1], a[0][1] = a[1][0], a[1][0] = x;
  x = a[0][2], a[0][2] = a[2][0], a[2][0] = x;
  x = a[2][1], a[2][1] = a[1][2], a[1][2] = x;
  return a;
}

/* a = u^T v */
INLINE dv3_t *dm3_vtv(double a[3][3], const double *u, const double *v)
{
  a[0][0] = u[0]*v[0];
  a[0][1] = u[0]*v[1];
  a[0][2] = u[0]*v[2];
  a[1][0] = u[1]*v[0];
  a[1][1] = u[1]*v[1];
  a[1][2] = u[1]*v[2];
  a[2][0] = u[2]*v[0];
  a[2][1] = u[2]*v[1];
  a[2][2] = u[2]*v[2];
  return a;
}

/* a += b */
INLINE dv3_t *dm3_inc(double a[3][3], double b[3][3])
{
  a[0][0] += b[0][0];
  a[0][1] += b[0][1];
  a[0][2] += b[0][2];
  a[1][0] += b[1][0];
  a[1][1] += b[1][1];
  a[1][2] += b[1][2];
  a[2][0] += b[2][0];
  a[2][1] += b[2][1];
  a[2][2] += b[2][2];
  return a;
}

/* a += b*s */
INLINE dv3_t *dm3_sinc(double a[3][3], double b[3][3], double s)
{
  a[0][0] += b[0][0]*s;
  a[0][1] += b[0][1]*s;
  a[0][2] += b[0][2]*s;
  a[1][0] += b[1][0]*s;
  a[1][1] += b[1][1]*s;
  a[1][2] += b[1][2]*s;
  a[2][0] += b[2][0]*s;
  a[2][1] += b[2][1]*s;
  a[2][2] += b[2][2]*s;
  return a;
}

/* c = a b */
INLINE dv3_t *dm3_mul(double c[3][3], double a[3][3], double b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
  return c;
}

/* c = a b^T */
INLINE dv3_t *dm3_mult(double c[3][3], double a[3][3], double b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = dv3_dot(a[i], b[j]);
  return c;
}

/* c = a^T b */
INLINE dv3_t *dm3_tmul(double c[3][3], double a[3][3], double b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[0][i]*b[0][j] + a[1][i]*b[1][j] + a[2][i]*b[2][j];
  return c;
}

/* c = a v */
INLINE double *dm3_mulvec(double *c, double a[3][3], const double *v)
{
  c[0] = a[0][0]*v[0] + a[0][1]*v[1] + a[0][2]*v[2];
  c[1] = a[1][0]*v[0] + a[1][1]*v[1] + a[1][2]*v[2];
  c[2] = a[2][0]*v[0] + a[2][1]*v[1] + a[2][2]*v[2];
  return c;
}

/* c = a^T v */
INLINE double *dm3_tmulvec(double *c, double a[3][3], const double *v)
{
  c[0] = a[0][0]*v[0] + a[1][0]*v[1] + a[2][0]*v[2];
  c[1] = a[0][1]*v[0] + a[1][1]*v[1] + a[2][1]*v[2];
  c[2] = a[0][2]*v[0] + a[1][2]*v[1] + a[2][2]*v[2];
  return c;
}

/* determinant of a 3x3 matrix */
INLINE double dm3_det(double a[3][3])
{
  return a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
      +  a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a[2][2])
      +  a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]);
}

/* inverse matrix b = a^(-1) */
INLINE dv3_t *dm3_inv(double b[3][3], double a[3][3])
{
  double d00, d01, d02, detm;
  d00 = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  d01 = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  d02 = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  detm = a[0][0]*d00 + a[0][1]*d01 + a[0][2]*d02;
  if (fabs(detm) < DBL_EPSILON) detm = (detm < 0) ? -DBL_EPSILON: DBL_EPSILON;
  b[0][0] = d00/detm;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detm;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/detm;
  b[1][0] = d01/detm;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/detm;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/detm;
  b[2][0] = d02/detm;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/detm;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detm;
  return b;
}

INLINE rv3_t *rm3_fromdm3(real a[3][3], double da[3][3])
{
  a[0][0] = (real) da[0][0];
  a[0][1] = (real) da[0][1];
  a[0][2] = (real) da[0][2];
  a[1][0] = (real) da[1][0];
  a[1][1] = (real) da[1][1];
  a[1][2] = (real) da[1][2];
  a[2][0] = (real) da[2][0];
  a[2][1] = (real) da[2][1];
  a[2][2] = (real) da[2][2];
  return a;
}

/* a = b */
INLINE rv3_t *rm3_copy(real a[3][3], real b[3][3])
{
  rv3_copy(a[0], b[0]);
  rv3_copy(a[1], b[1]);
  rv3_copy(a[2], b[2]);
  return a;
}

/* transpose */
INLINE rv3_t *rm3_trans(real a[3][3]) 
{
  real x;
  x = a[0][1], a[0][1] = a[1][0], a[1][0] = x;
  x = a[0][2], a[0][2] = a[2][0], a[2][0] = x;
  x = a[2][1], a[2][1] = a[1][2], a[1][2] = x;
  return a;
}

/* a = u^T v */
INLINE rv3_t *rm3_vtv(real a[3][3], const real *u, const real *v)
{
  a[0][0] = u[0]*v[0];
  a[0][1] = u[0]*v[1];
  a[0][2] = u[0]*v[2];
  a[1][0] = u[1]*v[0];
  a[1][1] = u[1]*v[1];
  a[1][2] = u[1]*v[2];
  a[2][0] = u[2]*v[0];
  a[2][1] = u[2]*v[1];
  a[2][2] = u[2]*v[2];
  return a;
}

/* a += b */
INLINE rv3_t *rm3_inc(real a[3][3], real b[3][3])
{
  a[0][0] += b[0][0];
  a[0][1] += b[0][1];
  a[0][2] += b[0][2];
  a[1][0] += b[1][0];
  a[1][1] += b[1][1];
  a[1][2] += b[1][2];
  a[2][0] += b[2][0];
  a[2][1] += b[2][1];
  a[2][2] += b[2][2];
  return a;
}

/* a += b*s */
INLINE rv3_t *rm3_sinc(real a[3][3], real b[3][3], real s)
{
  a[0][0] += b[0][0]*s;
  a[0][1] += b[0][1]*s;
  a[0][2] += b[0][2]*s;
  a[1][0] += b[1][0]*s;
  a[1][1] += b[1][1]*s;
  a[1][2] += b[1][2]*s;
  a[2][0] += b[2][0]*s;
  a[2][1] += b[2][1]*s;
  a[2][2] += b[2][2]*s;
  return a;
}

/* c = a b */
INLINE rv3_t *rm3_mul(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
  return c;
}

/* c = a b^T */
INLINE rv3_t *rm3_mult(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = rv3_dot(a[i], b[j]);
  return c;
}

/* c = a^T b */
INLINE rv3_t *rm3_tmul(real c[3][3], real a[3][3], real b[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      c[i][j] = a[0][i]*b[0][j] + a[1][i]*b[1][j] + a[2][i]*b[2][j];
  return c;
}

/* c = a v */
INLINE real *rm3_mulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0] + a[0][1]*v[1] + a[0][2]*v[2];
  c[1] = a[1][0]*v[0] + a[1][1]*v[1] + a[1][2]*v[2];
  c[2] = a[2][0]*v[0] + a[2][1]*v[1] + a[2][2]*v[2];
  return c;
}

/* c = a^T v */
INLINE real *rm3_tmulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0] + a[1][0]*v[1] + a[2][0]*v[2];
  c[1] = a[0][1]*v[0] + a[1][1]*v[1] + a[2][1]*v[2];
  c[2] = a[0][2]*v[0] + a[1][2]*v[1] + a[2][2]*v[2];
  return c;
}

/* determinant of a 3x3 matrix */
INLINE real rm3_det(real a[3][3])
{
  return a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
      +  a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a[2][2])
      +  a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]);
}

/* inverse matrix b = a^(-1) */
INLINE rv3_t *rm3_inv(real b[3][3], real a[3][3])
{
  real d00, d01, d02, detm;
  d00 = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  d01 = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  d02 = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  detm = a[0][0]*d00 + a[0][1]*d01 + a[0][2]*d02;
  if (fabs(detm) < 1e-30) detm = (detm < 0) ? -1e-30f: 1e-30f;
  b[0][0] = d00/detm;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detm;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/detm;
  b[1][0] = d01/detm;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/detm;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/detm;
  b[2][0] = d02/detm;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/detm;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detm;
  return b;
}

/* compute eigenvalues of a 3x3 matrix
 * solving a cubic equation */
INLINE double *dm3_eigval(double v[3], double a[3][3])
{
  double m, p, q, pr, pr3, a00, a11, a22;

  m = (a[0][0] + a[1][1] + a[2][2])/3;
  a00 = a[0][0] - m;
  a11 = a[1][1] - m;
  a22 = a[2][2] - m;
  q = ( a00 * (a11*a22 - a[1][2]*a[2][1])
      + a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a22)
      + a[0][2] * (a[1][0]*a[2][1] - a11*a[2][0]) ) / 2.0;
  p = (a00*a00 + a11*a11 + a22*a22) / 6.0
    + (a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]) / 3.0;
  /* solve x^3 - 3 p x  - 2 q = 0 */
  pr = sqrt(p);
  pr3 = p*pr;
  if (pr3 <= fabs(q)) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = m + pr;
      v[2] = m - 2.0 * pr;
    } else { /* phi = 0 */
      v[0] = m + 2.0 * pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    double phi = acos(q/pr3)/3.0; /* 0 < phi < pi/3 */
    v[0] = m + 2.0 * pr * cos(phi);  /* largest */
    v[1] = m + 2.0 * pr * cos(phi - 2*M_PI/3); /* second largest */
    v[2] = m + 2.0 * pr * cos(phi + 2*M_PI/3); /* smallest */
#ifdef RV3_DEBUG
    { int i; double vi, y[3], dy[3];
      for (i = 0; i < 3; i++) { vi = v[i] - m; y[i] = vi*(vi*vi - 3*p) - 2 * q; dy[i] = 3*(vi*vi - p); }
      dv3_print(v,  "roots   ", "%26.14e", 1);
      dv3_print(y,  "residues", "%26.14e", 1);
      dv3_print(dy, "slope   ", "%26.14e", 1);
      printf("m %20.14f, q/pr3 %22.14e, pr %20.14f, phi %20.14f, %20.14f deg\n", m, q/pr3, pr, phi/M_PI*180, 120 - phi/M_PI*180);
    }
#endif
  }
  return v;
}

/* sort s to descending order, order u and v correspondingly */
INLINE void dv3_sort3(double s[3], double (*u)[3], double (*v)[3])
{
  double tmp;

  if (s[2] > s[1]) { 
    tmp = s[1]; s[1] = s[2]; s[2] = tmp;
    if (u) dv3_swap(u[1], u[2]);
    if (v) dv3_swap(v[1], v[2]);
  }
  if (s[1] > s[0]) {
    tmp = s[0]; s[0] = s[1]; s[1] = tmp;
    if (u) dv3_swap(u[0], u[1]);
    if (v) dv3_swap(v[0], v[1]);
  }
  if (s[2] > s[1]) {
    tmp = s[1]; s[1] = s[2]; s[2] = tmp;
    if (u) dv3_swap(u[1], u[2]);
    if (v) dv3_swap(v[1], v[2]);
  }
}

/* return the pivot row for column c, row index starts from r0 */
INLINE int dm3_pivot_(double m[3][3], int r0, int c, double *max)
{ 
  int i, r = r0; 
  double tmp;

  for (*max = fabs(m[r = r0][c]), i = r0 + 1; i < 3; i++) 
    if ((tmp = fabs(m[i][c])) > *max) r = i, *max = tmp;
  return r;
}
 
/* solve matrix equation a x = 0, matrix 'a' is destroyed
 * solutions are saved as *row* vectors in 'x'
 * return the number of solutions */
INLINE int dm3_solvezero(double a[3][3], double (*x)[3], double tol)
{
  double max;
  int i, j, k, ns = 0;

  k = dm3_pivot_(a, 0, 0, &max); /* pivot for column 0 */
  if (max <= tol) { /* found the first eigenvector */
    dv3_make(x[ns++], 1, 0, 0);
    k = dm3_pivot_(a, 0, 1, &max); /* pivot for column 1 */
    if (max <= tol) {
      dv3_make(x[ns++], 0, 1, 0);
      k = dm3_pivot_(a, 0, 2, &max);
      if (max <= tol) dv3_make(x[ns++], 0, 0, 1);
    } else {
      if (k != 0) dv3_swap(a[0], a[k]);
      a[0][2] /= a[0][1]; /* normalize row 0, a[0][1] = 1; */
      for (i = 1; i < 3; i++) a[i][2] -= a[i][1]*a[0][2];
      dm3_pivot_(a, 1, 2, &max);
      if (max <= tol) dv3_makenorm(x[ns++], 0, a[0][2], -1);
    }
  } else {
    if (k != 0) dv3_swap(a[0], a[k]);
    a[0][1] /= a[0][0]; a[0][2] /= a[0][0]; /* normalize row 0, a[0][0] = 1 */
    for (i = 1; i < 3; i++)
      for (j = 1; j < 3; j++) /* a[i][0] = 0 now */
        a[i][j] -= a[i][0]*a[0][j];
    k = dm3_pivot_(a, 1, 1, &max); /* pivot for column 1 */
    if (max <= tol) { /* column 1 is empty */
      dv3_makenorm(x[ns++], -a[0][1], 1, 0);
      k = dm3_pivot_(a, 1, 2, &max);
      if (max <= tol) /* column 2 is empty too */
        dv3_makenorm(x[ns++], -a[0][2], 0, 1);
    } else {
      if (k != 1) dv3_swap(a[1], a[k]);
      a[1][2] /= a[1][1]; /* normalize row 1, a[1][1] = 1 */
      a[2][2] -= a[2][1]*a[1][2];
      if (fabs(a[2][2]) > tol) {
#ifdef RV3_DEBUG
        printf("a22 %g vs tol %g\n", a[2][2], tol);
#endif
        return 0; /* no solutions */
      }
      a[0][2] -= a[0][1]*a[1][2];
      dv3_makenorm(x[ns++], -a[0][2], -a[1][2], 1);
    }
  }
  return ns;
}

/* given an eigenvalue, return the corresponding eigenvectors */
INLINE int dm3_eigvecs(double (*vecs)[3], double mat[3][3], double val, double tol)
{
  double m[3][3];

  dm3_copy(m, mat); /* make a matrix */
  m[0][0] -= val; m[1][1] -= val; m[2][2] -= val;
  return dm3_solvezero(m, vecs, tol);
}

/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs' */
INLINE dv3_t *dm3_eigsys(double v[3], double vecs[3][3], double mat[3][3], int nt)
{
  double vs[5][3], sq, tol, nv; /* for safty, vs needs 5 rows */
  int n = 0, nn, i;
 
  dm3_eigval(v, mat);
  for (sq = 0, i = 0; i < 3; i++) sq += dv3_sqr(mat[i]);
  /* errors of the eigenvalues from the cubic equation can reach sqrt(eps)
   * use a large tolerance */
  tol = 10.0 * sqrt(sq * DBL_EPSILON);

  for (nn = i = 0; i < 3; i++) {
    n = dm3_eigvecs(vs+nn, mat, v[nn], tol);
    if (n == 0) goto ERR;
    if ((nn += n) >= 3) break;
  }
  
  /* NOTE: make sure eigenvectors are orthogonal */
  dv3_normalize( dv3_cross(vs[2], vs[0], vs[1]) );
  dv3_normalize( dv3_cross(vs[1], vs[2], vs[0]) );

  dm3_copy(vecs, vs);
  for (i = 0; i < 3; i++) {
    nv = dv3_dot(dm3_mulvec(vs[i], mat, vecs[i]), vecs[i]);
    if (fabs(nv - v[i]) > tol) {
      fprintf(stderr, "corrupted eigenvalue i %d, %g vs. %g\n", i, nv, v[i]);
      goto ERR;
    }
#ifdef RV3_DEBUG
    printf("Eigenvalue: %22.14f vs %22.14f (corrected)\n", v[i], nv);
    dv3_print(vecs[i], "eigenvector i", "%20.12e", 1);
#endif
    v[i] = nv;
  }
#ifdef RV3_DEBUG
  printf("det(V) = %g\n", dm3_det(vecs));
#endif
  dv3_sort3(v, vecs, NULL);

  if (nt) return vecs; else return dm3_trans(vecs);
ERR:
  printf("fatal: bad eigenvalues, n %d, nn %d\n", n, nn);
  dm3_print(mat, "matrix", "%24.16e", 1);
  dv3_print(v, "eigenvalues", "%24.16e", 1);
  exit(1);
  return NULL;
}

/* SVD decomposition of a 3x3 matrix A = U S V^T */
INLINE void dm3_svd(double a[3][3], double u[3][3], double s[3], double v[3][3])
{
  int i, rank;
  double ata[3][3], us[3][3];

  /* 1. compute A^T A and its eigenvectors, which is V */
  dm3_tmul(ata, a, a);
  dm3_eigsys(s, v, ata, 1);
#ifdef RV3_DEBUG
  dv3_print(s, "S^2 ", "%22.14e", 1);
  dm3_print(ata, "A^T A ",  "%20.14f", 1);
  dm3_print(v, "V^T ",  "%20.14f", 1);
#endif

  /* 2. U^T = S^{-1} V^T A^T, and each row of U^T is an eigenvector
   * since eigenvectors are to be normalized, S^{-1} is unnecessary */
  if (s[0] <= 0.0) {
    rank = 0;
    dm3_copy(u, v);
  } else {
    double tol = 10. * sqrt(DBL_EPSILON);
    /* the test i = 1 + (s[1] > s[0]*tol) + (s[2] > s[0]*tol); */
    dm3_mult(u, v, a);
    for (i = 0; i < 3; i++) {
      dv3_copy(us[i], u[i]); /* save a copy of V^T A^T before normalizing it */
      s[i] = dv3_norm(u[i]);
      if (s[i] > 0.0) dv3_smul(u[i], 1.0/s[i]);
    }
    rank = 1;
    rank += (fabs(dv3_dot(u[0], u[1])) < tol && s[1] > tol);
    rank += (fabs(dv3_dot(u[0], u[2])) < tol && fabs(dv3_dot(u[1], u[2])) < tol && s[2] > tol);
#ifdef RV3_DEBUG
    dm3_print(u, "U^T ", "%22.14e", 1);
    dm3_print(us, "Us^T ", "%22.14e", 1);
    dv3_print(s, "S ", "%22.14e", 1);
    dm3_print(a, "A ",  "%20.14f", 1);
    printf("rank = %d, tol %g, det %g, u0.u0 %g, u0.u1 %g, u0.u2 %g, u1.u1 %g, u1.u2 %g, u2.u2 %g\n\n\n", rank, tol, dm3_det(u), 
        dv3_sqr(u[0]), dv3_dot(u[0], u[1]), dv3_dot(u[0], u[2]), dv3_sqr(u[1]), dv3_dot(u[1], u[2]), dv3_sqr(u[2]));
#endif
    if (rank <= 2) {
      if (rank == 1) {
        double z[3] = {0, 0, 0}, w, tmp;
        w = fabs(u[0][i = 0]);
        if ((tmp = fabs(u[0][1])) < w) w = tmp, i = 1;
        if ((tmp = fabs(u[0][2])) < w) i = 2;
        z[i] = 1.0f; /* select the smallest element in u[0] as z */
        dv3_normalize( dv3_cross(u[1], z, u[0]) );
        s[1] = dv3_dot(u[1], us[1]); /* S = U^T (V^T A^T)^T is more accurate than sqrt(A^T A) */
        if (s[1] < 0) { s[1] = -s[1]; dv3_neg(u[1]); } /* make sure s[1] > 0 */
      }
      dv3_normalize( dv3_cross(u[2], u[0], u[1]) );
      s[2] = dv3_dot(u[2], us[2]);
      if (s[2] < 0) { s[2] = -s[2]; dv3_neg(u[2]); }
    }
    dv3_sort3(s, u, v);
#ifdef RV3_DEBUG
    printf("det(U) %g, det(V) %g\n", dm3_det(u), dm3_det(v));
#endif
  }
  dm3_trans(v);
  dm3_trans(u);
}

/* eigenvalues of a 3x3 matrix 
 * internal calculation are carried out in double precision */
INLINE real *rm3_eigval(real v[3], real a[3][3])
{
  if (sizeof(real) == sizeof(double)) {
    return (real *) dm3_eigval((double *) v, (dv3_t *) a);
  } else { /* the routine require high precision, double is safer */
    double da[3][3], dv[3];

    dm3_fromrm3(da, a);
    dm3_eigval(dv, da);
    return rv3_fromdv3(v, dv);
  } 
}

/* solve A x = 0 */
INLINE int rm3_solvezero(real a[3][3], real (*x)[3], real tol)
{
  if (sizeof(real) == sizeof(double)) {
    return dm3_solvezero((dv3_t *) a, (dv3_t *) x, tol);
  } else {
    double da[3][3], dx[3][3];
    int n, i;
    dm3_fromrm3(da, a);
    n = dm3_solvezero(da, dx, tol);
    for (i = 0; i < n; i++)
      rv3_fromdv3(x[i], dx[i]);
    return n;
  }
}

/* given an eigenvalue, return the corresponding eigenvectors */
INLINE int rm3_eigvecs(real (*vecs)[3], real mat[3][3], real val, real tol)
{
  real m[3][3];

  rm3_copy(m, mat); /* make a matrix */
  m[0][0] -= val; m[1][1] -= val; m[2][2] -= val;
  return rm3_solvezero(m, vecs, tol);
}

/* compute eigenvectors for the eigenvalues
 * ideally, eigenvalues should be sorted in magnitude-descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
INLINE rv3_t *rm3_eigsys(real v[3], real vecs[3][3], real mat[3][3], int nt)
{
  if (sizeof(real) == sizeof(double)) {
    return (rv3_t *) dm3_eigsys((double *) v, (dv3_t *) vecs, (dv3_t *) mat, nt);
  } else {
    double dvecs[3][3], dmat[3][3], dv[3];

    dm3_fromrm3(dmat, mat);
    dm3_eigsys(dv, dvecs, dmat, nt);
    rv3_fromdv3(v, dv);
    return rm3_fromdm3(vecs, dvecs);
  } 
}

/* SVD decomposition of a 3x3 matrix a = u s v^T */
INLINE void rm3_svd(real a[3][3], real u[3][3], real s[3], real v[3][3])
{
  if (sizeof(real) == sizeof(double)) {
    dm3_svd((dv3_t *) a, (dv3_t *) u, (double *) s, (dv3_t *) v);
  } else {
    double da[3][3], du[3][3], ds[3], dv[3][3];
    
    dm3_fromrm3(da, a);
    dm3_svd(da, du, ds, dv);
    rm3_fromdm3(u, du);
    rm3_fromdm3(v, dv);
    rv3_fromdv3(s, ds);
  }
}

/* return 0 rotation matrix around v for ang */
INLINE rv3_t *rm3_mkrot(real m[3][3], const real *v, real ang)
{
  real c = (real) cos(ang), s = (real) sin(ang), nc, n[3];
  
  rv3_copy(n, v);
  rv3_normalize(n);
  nc = 1 - c;
  m[0][0] = n[0]*n[0]*nc + c;
  m[0][1] = n[0]*n[1]*nc - n[2]*s;
  m[0][2] = n[0]*n[2]*nc + n[1]*s;
  m[1][0] = n[1]*n[0]*nc + n[2]*s;
  m[1][1] = n[1]*n[1]*nc + c;
  m[1][2] = n[1]*n[2]*nc - n[0]*s;
  m[2][0] = n[2]*n[0]*nc - n[1]*s;
  m[2][1] = n[2]*n[1]*nc + n[0]*s;
  m[2][2] = n[2]*n[2]*nc + c;
  return m;
}

/* rotate v0 around u by ang, save result to v1 */
INLINE real *rv3_rot(real *v1, const real *v0, const real *u, real ang)
{
  real m[3][3];

  rm3_mkrot(m, u, ang);
  rm3_mulvec(v1, m, v0);
  return v1;
}

/* uniformly distributed random vector [a, a + b) */
#define rv3_rnd0() rv3_rnd(v, 0, 1)
INLINE real *rv3_rnd(rv3_t v, real a, real b)
{
  v[0] = (real) (a + b * rnd0());
  v[1] = (real) (a + b * rnd0());
  v[2] = (real) (a + b * rnd0());
  return v;
}

/* normally distributed random vector */
#define rv3_grand0(v) rv3_grand(v, 0, 1)
INLINE real *rv3_grand(rv3_t v, real c, real r)
{
  v[0] = (real) (c + r * grand0());
  v[1] = (real) (c + r * grand0());
  v[2] = (real) (c + r * grand0());
  return v;
}

/* generate a random orthonormal (unitary) 3x3 matrix */
INLINE rv3_t *rm3_rnduni(real a[3][3])
{
  real dot;

  rv3_rnd(a[0], -.5f, 1.f);
  rv3_normalize(a[0]);

  rv3_rnd(a[1], -.5f, 1.f);
  /* normalize a[1] against a[0] */
  dot = rv3_dot(a[0], a[1]);
  rv3_sinc(a[1], a[0], -dot);
  rv3_normalize(a[1]);

  rv3_cross(a[2], a[0], a[1]);
  return a;
}

#endif /* RV3_H__ */

