#ifndef INLINE
#define INLINE __inline static
#endif
#define RESTRICT __restrict 
#include "def.h"
#ifndef RV3_H__
#define RV3_H__

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

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */

INLINE real *rv3_make(real *x, real a, real b, real c)
  { x[0] = a; x[1] = b; x[2] = c; return x; }
INLINE real *rv3_zero(real *x) { return rv3_make(x, 0, 0, 0); }
INLINE real *rv3_copy(real *x, const real *src)
  { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; return x; }
/* use macro to avoid const qualifier of src */
#define rv3_ncopy(x, src, n) memcpy(x, src, n*sizeof(x[0]))

INLINE real rv3_sqr (const real *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
INLINE real rv3_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

/* if x == y, try to use sqr */
INLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

INLINE real *rv3_cross(real *RESTRICT z, const real *x, const real *y)
{
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
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
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  x[2] += s*dx[2];
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
  if (r > 0.0) rv3_smul(x, 1.f/r);
  return x;
}

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
  sgn = ((vol > 0.0f) ? 1.0f : (-1.0f));
  phi *= sgn;
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
      real vgi[3], vgj[3], vgk[3], vgl[3];
      real uvec[3], vvec[3], svec[3], p, q;
      real gi2, gj2, gk2, gl2, g2all, invg2;
      unsigned doi, doj, dok, dol;

      doi = (flags & DIH_I);
      doj = (flags & DIH_J);
      dok = (flags & DIH_K);
      dol = (flags & DIH_L);

      scl = nxkj/m2;
      rv3_smul2(vgi, m, scl);
      scl = -nxkj/n2;
      rv3_smul2(vgl, n, scl);

      p = rv3_dot(xij, xkj);
      p /= nxkj2;
      rv3_smul2(uvec, vgi, p);
      q = rv3_dot(xkl, xkj);
      q /= nxkj2;
      rv3_smul2(vvec, vgl, q);
      rv3_diff(svec, uvec, vvec);

      rv3_diff(vgj, svec, vgi);
      rv3_nadd(vgk, vgl, svec);

      rv3_copy(dih->g[0], vgi);
      rv3_copy(dih->g[1], vgj);
      rv3_copy(dih->g[2], vgk);
      rv3_copy(dih->g[3], vgl);

      gi2 = rv3_sqr(vgi);
      gj2 = rv3_sqr(vgj);
      gk2 = rv3_sqr(vgk);
      gl2 = rv3_sqr(vgl);
      g2all = 0.0f;
      if (doi) g2all += gi2;
      if (doj) g2all += gj2;
      if (dok) g2all += gk2;
      if (dol) g2all += gl2;
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

        gjxij = rv3_dot(vgj, xij);
        gjxkl = rv3_dot(vgj, xkl);
        gjmvv = rv3_dot(vgj, mvv);
        gjnvv = rv3_dot(vgj, nvv);
        gkxij = rv3_dot(vgk, xij);
        gkxkl = rv3_dot(vgk, xkl);
        gkmvv = rv3_dot(vgk, mvv);
        gknvv = rv3_dot(vgk, nvv);

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
      c[i][j] = a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
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

/* c = a v */
INLINE real *rm3_mulvec(real *c, real a[3][3], const real *v)
{
  c[0] = a[0][0]*v[0] + a[0][1]*v[1] + a[0][2]*v[2];
  c[1] = a[1][0]*v[0] + a[1][1]*v[1] + a[1][2]*v[2];
  c[2] = a[2][0]*v[0] + a[2][1]*v[1] + a[2][2]*v[2];
  return c;
}

/* c = a^T v */
INLINE real *rm3_multvec(real *c, real a[3][3], const real *v)
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

/* eigenvalues of a 3x3 matrix */
INLINE real *rm3_eigval(real v[3], real a[3][3])
{
  real m, p, q, cphi, sphi, pr, pr3;

  m = (a[0][0]+a[1][1]+a[2][2])/3.f;
  a[0][0] -= m;
  a[1][1] -= m;
  a[2][2] -= m;
  q = .5f * rm3_det(a);
  p = ((a[0][0]*a[0][0] + a[1][1]*a[1][1] + a[2][2]*a[2][2]) +
   2.f*(a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]))/6.f;
  pr = (real) sqrt(p);
  pr3 = p*pr;
  if (pr3 <= fabs(q)) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = m + pr;
      v[2] = m - 2.f*pr;
    } else { /* phi = 0 */
      v[0] = m + 2.f*pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    double phi = acos(q/pr3)/3.f; /* 0 < phi < pi/3 */
    cphi = (real)cos(phi);
    sphi = (real)(sin(phi)*1.7320508075688772);
    v[0] = m + 2.f*pr*cphi;  /* cos(phi), largest */
    v[1] = m - pr*(cphi-sphi); /* cos(phi-2*pi/3), second largest */
    v[2] = m - pr*(cphi+sphi); /* cos(phi+2*pi/3), smallest */
  }
  a[0][0] += m;
  a[1][1] += m;
  a[2][2] += m;  
  return v;
}

/* given matrix a and eigenvalue lm, return eigenvector */
INLINE real *rm3_eigvec(real vec[3], real m[3][3], real val)
{
  double a = m[0][0]-val, b = m[1][1]-val, c, d = m[0][1], e = m[0][2], f = m[1][2];
  double detm, tol = 1e-12;

  vec[2] = 1.f;
  if (fabs(detm = a*b - d*d) > tol) { /* use row 0 and 1 */
    rv3_make(vec, (real)((d*f-b*e)/detm), (real)((e*d-a*f)/detm), 1);
    return rv3_normalize(vec);
  }
  c = m[2][2] - val;
  if (fabs(detm = a*f - e*d) > tol) { /* row 1 and 2 */
    rv3_make(vec, (real)((d*c-e*f)/detm), (real)((e*e-a*c)/detm), 1);
    return rv3_normalize(vec);
  }
  if ((detm = sqrt(a*a+d*d)) > tol) { /* three-row-degenerate */
    rv3_make(vec, (real)(d/detm), (real)(-a/detm), 0);
  } else {
    rv3_make(vec, 1, 0, 0);
  }
  return vec;
}

/* compute eigenvectors for the eigenvalues */
INLINE rv3_t *rm3_eigvecs(real vecs[3][3], real mat[3][3], real v[3], int t)
{
  const double tol = 1e-12;
  double v0 = fabs(v[0]), v1 = fabs(v[1]), v2 = fabs(v[2]);
  
  rm3_eigvec(vecs[0], mat, v[0]);
  if ( fabs(v[0] - v[1]) > tol*(v0 + v1) ) {
    rm3_eigvec(vecs[1], mat, v[1]);
    rv3_cross(vecs[2], vecs[0], vecs[1]);
  } else if ( fabs(v[2] - v[1]) > tol*(v1 + v2) ) {
    rm3_eigvec(vecs[2], mat, v[2]);
    rv3_cross(vecs[1], vecs[2], vecs[0]);
  } else {
    rv3_make(vecs[1], 0, 1, 0);
    rv3_make(vecs[2], 0, 0, 1);
  }
  /* transpose the matrix */
  if (t) return vecs;
  else return rm3_trans(vecs);
}

/* SVD decomposition of a 3x3 matrix a = u s v */
INLINE int rm3_svd(real a[3][3], real u[3][3], real s[3], real v[3][3])
{
  int i, j;
  real ata[3][3], z[3];
  const double tol = 1e-12;

  /* 1. compute A^T A and its eigenvectors */
  for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) {
      ata[i][j] = a[0][i]*a[0][j] + a[1][i]*a[1][j] + a[2][i]*a[2][j];
      if (i != j) ata[j][i] = ata[i][j];
    }
  rm3_eigval(s, ata);
  for (i = 0; i < 3; i++) if (s[i] < 0) s[i] = 0;
  rm3_eigvecs(v, ata, s, 1); /* get V^T */

  /* 2. U = A V S^-1, or U^T = S^{-1}T V^T A^T */
  j = (s[0] > tol) + (s[1] > tol) + (s[2] > tol);
  if (j >= 2) {
    rm3_mult(u, v, a);
    if (j == 2) rv3_cross(u[2], u[0], u[1]); /* fix the last */
  } else if (j == 1) {
    rm3_multvec(u[0], a, v[0]);
    rv3_zero(z);
    /* choose z[i] such that z X u[0] != 0 */
    i = (u[0][0]*u[0][0] < u[0][1]*u[0][1]) ? 0 : 1;
    rv3_cross(u[1], z, u[0]);
    rv3_cross(u[2], u[0], u[1]);
  } else { /* use u */
    for (i = 0; i < 3; i++) rv3_copy(u[i], v[i]);
  }
  for (i = 0; i < 3; i++) rv3_normalize(u[i]);
  for (i = 0; i < 3; i++) s[i] = (real) sqrt(s[i]);
  rm3_trans(u);
  rm3_trans(v);
  return 0;
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

#define rv3_print(r, nm, fmt, nl) rv3_fprint(stdout, r, nm, fmt, nl)
INLINE void rv3_fprint(FILE *fp, const real *r, const char *nm,
    const char *fmt, int nl)
{
  int i;
  if (nm) fprintf(fp, "%s: ", nm);
  for (i = 0; i < 3; i++)
    fprintf(fp, fmt, r[i], nl);
  fprintf(fp, "%c", (nl ? '\n' : ';'));
}

#define rm3_print(r, nm, fmt, nl) rm3_fprint(stdout, r, nm, fmt, nl)
INLINE void rm3_fprint(FILE *fp, real r[3][3], const char *nm,
    const char *fmt, int nl)
{
  int i, j;
  if (nm) fprintf(fp, "%s:%c", nm, (nl ? '\n' : ' '));
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      fprintf(fp, fmt, r[i][j], nl);
    }
    fprintf(fp, "%s", (nl ? "\n" : "; "));
  }
}

#endif /* RV3_H__ */

