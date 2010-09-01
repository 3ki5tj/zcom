#ifndef DIHCALC_C__
#define DIHCALC_C__

#include <stdio.h>
#include <math.h>
#include "dihcalc.h"

/* Calculates the dihedral angle, its gradient and the divegence
   of a conjugate field (to the gradient).

   For simplest calculation, dc can be NULL and flags=0,
   in this case only the dihedral angle is computed.

   The gradient and divergent are optional, depending on the flags,
   DIHCALC_GRAD and DIHCALC_DIV.
   To calculate the divergence, the force must be calculated.

   If you want to use a special routine for treating
   periodic boundary condition, assign a function pointer
   to pbc_rv3_diff before calling; and also supply additional
   information by the `pbc' entry.
   Otherwise, pbc_rv3_diff *must* be set to NULL.

   The procedures of calculating the angle and force
   are similar to that in GROMACS.

   We also introduce the conjugate field v
     v = grad(phi) / [grad(phi).grad(phi)],
   such that v.grad(phi) = 1.0.
   This field is not explicitly computed (since it's just a simple
   scaling of the gradient).
   The denominator is saved to grad2;

   You include only a few components of the gradient in calculating
   the conjugate field by passing `flags' a combination of
   DIHCALC_I, DIHCALC_J, DIHCALC_K, and DIHCALC_L.

   We also calculate the divergence of v, or div.v,
   It turns out the divergence of the gradient itself is always zero,
     div.grad(phi) = 0.
   Thus, we have
     div.v = -2 [ grad(phi).(grad grad(phi)).grad(phi) ] /[grad(phi).grad(phi)]^2.
   Components in grad(phi).(grad grad(phi)).grad(phi) involve the
   four gradient terms: d4ij, d4ik, ..., d4kl.

   Note, no matter what combination of DIHCALC_I, DIHCALC_J, DIHCALC_K, and DIHCALC_L
   is used, *all* moments d4ij, d4ik, ... , d4kl are calculated for safety.
   But only the involved ones are used to combined to produce the divergence.
 */
real rv3_calcdih(dihcalc_t *dc,
    const real *xi, const real *xj, const real *xk, const real *xl,
    unsigned int flags)
{
  real dot, scl, tol, vol, phi, sign, cosphi;
  real nxkj, nxkj2, m2, n2;
  real xij[3], xkj[3], xkl[3];
  real m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */
  const real cosmax = (real)(1.0-6e-8);

  if (dc != NULL && dc->pbc_rv3_diff != NULL) { /* handle pbc */
    dc->t1 = (*dc->pbc_rv3_diff)(dc->pbc, xi, xj, xij);
    dc->t2 = (*dc->pbc_rv3_diff)(dc->pbc, xk, xj, xkj);
    dc->t3 = (*dc->pbc_rv3_diff)(dc->pbc, xk, xl, xkl);
  } else {
    rv3_diff(xij, xi, xj);
    rv3_diff(xkj, xk, xj);
    rv3_diff(xkl, xk, xl);
  }
  nxkj2 = rv3_sqr(xkj);
  nxkj = (real)sqrt(nxkj2);
  tol = nxkj2*6e-8f;

  rv3_cross(m, xij, xkj);
  m2 = rv3_sqr(m);
  rv3_cross(n, xkj, xkl);
  n2 = rv3_sqr(n);
  if (m2 > tol && n2 > tol) {
    scl = (real)sqrt(m2*n2);
    dot = rv3_dot(m, n);
    cosphi = dot/scl;
    if (cosphi > cosmax) cosphi = cosmax; else if (cosphi < -cosmax) cosphi = -cosmax;
  } else {
    cosphi = cosmax;
  }
  phi = (real)acos(cosphi);
  vol = rv3_dot(n, xij);
  sign = ((vol > 0.0f)?1.0f:(-1.0f));
  phi *= sign;
  if (dc != NULL) {
    dc->phi = phi;
    dc->sign = sign;
    dc->cosphi = cosphi;
    dc->flags = flags;
  }

  /* optionally calculate the gradient */
  if (dc != NULL && (flags&(DIHCALC_GRAD|DIHCALC_DIV)) ) { /* do gradient if divergence is required */
    /* clear divergence */
    dc->div = dc->d4ij = dc->d4ik = dc->d4jj = dc->d4jk = dc->d4jl = dc->d4kk = dc->d4kl = 0.0f;

    /* calculate the gradient of the force
     * the direction of which increase the dihedral by one unit.
     */

    if (m2 > tol && n2 > tol) {
      real gi[3], gj[3], gk[3], gl[3];
      real uvec[3], vvec[3], svec[3], p, q;
      real gi2, gj2, gk2, gl2, g2all, invg2;
      unsigned doi, doj, dok, dol;

      doi = (flags&DIHCALC_I);
      doj = (flags&DIHCALC_J);
      dok = (flags&DIHCALC_K);
      dol = (flags&DIHCALC_L);

      scl = nxkj/m2;
      rv3_smul2(gi, m, scl);
      scl = -nxkj/n2;
      rv3_smul2(gl, n, scl);

      p = rv3_dot(xij, xkj);
      p /= nxkj2;
      rv3_smul2(uvec, gi, p);
      q = rv3_dot(xkl, xkj);
      q /= nxkj2;
      rv3_smul2(vvec, gl, q);
      rv3_diff(svec, uvec, vvec);

      rv3_diff(gj, svec, gi);
      rv3_nsum2(gk, gl, svec);

      rv3_copy(dc->grad[0], gi);
      rv3_copy(dc->grad[1], gj);
      rv3_copy(dc->grad[2], gk);
      rv3_copy(dc->grad[3], gl);

      gi2 = rv3_sqr(gi);
      gj2 = rv3_sqr(gj);
      gk2 = rv3_sqr(gk);
      gl2 = rv3_sqr(gl);
      g2all = 0.0f;
      if (doi) g2all += gi2;
      if (doj) g2all += gj2;
      if (dok) g2all += gk2;
      if (dol) g2all += gl2;
      dc->grad2 = g2all;
      invg2 = 1.0f/g2all;

      if (flags&DIHCALC_DIV) {
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

        gjxij = rv3_dot(gj, xij);
        gjxkl = rv3_dot(gj, xkl);
        gjmvv = rv3_dot(gj, mvv);
        gjnvv = rv3_dot(gj, nvv);
        gkxij = rv3_dot(gk, xij);
        gkxkl = rv3_dot(gk, xkl);
        gkmvv = rv3_dot(gk, mvv);
        gknvv = rv3_dot(gk, nvv);

        tmp1 = nxkj2*sinmn;
        tmp2 = tmp1/m2;
        dc->d4ij = kikl*tmp2;
        dc->d4ik = ijlj*tmp2;
        tmp2 = tmp1/n2;
        dc->d4jl = kikl*tmp2;
        dc->d4kl = ijlj*tmp2;

        dc->d4jj = -(gjxij*gjmvv+gjxkl*gjnvv)/nxkj
                +2.0f*(kivkj*gjmvv-klvkj*gjnvv)*(-kikl*sinmn);

        dc->d4jk = (gjxij*gkmvv+gjxkl*gknvv)/nxkj
              +(-(gjmvv*ljvkj+gkmvv*klvkj)*(ijvkj*kivkj)
                +(gjnvv*ijvkj+gknvv*kivkj)*(ljvkj*klvkj) )*sinmn;

        dc->d4kk = -(gkxkl*gknvv+gkxij*gkmvv)/nxkj
                +2.0f*(ljvkj*gknvv-ijvkj*gkmvv)*(ijlj*sinmn);

        /* summarize */
        if ((flags&DIHCALC_FOUR) == DIHCALC_FOUR) {
          tmp1 = dc->d4jj + dc->d4kk;
          tmp2 = dc->d4ij + dc->d4ik+dc->d4jk+dc->d4jl+dc->d4kl;
        } else {
          tmp1 = tmp2 = 0.0f;
          if (doj) { tmp1 += dc->d4jj; }
          if (dok) { tmp1 += dc->d4kk; }
          if (doi && doj) tmp2 += dc->d4ij;
          if (doi && dok) tmp2 += dc->d4ik;
          if (doj && dok) tmp2 += dc->d4jk;
          if (doj && dol) tmp2 += dc->d4jl;
          if (dok && dol) tmp2 += dc->d4kl;
        }
        dc->div = -2.0f*(tmp1+2.0f*tmp2)*(invg2*invg2);
      } /* do divengence */

    } else { /* clear the gradients */
      int j;
      for (j = 0; j < 4; j++)
        rv3_zero(dc->grad[j]);
    }
  }

  return phi;
}

#endif

