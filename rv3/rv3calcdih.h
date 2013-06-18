#ifndef RV3CALCDIH_H__
#define RV3CALCDIH_H__
/* heavy-weight dihedral calculation */

/* structure for dihedral calculation */
typedef struct {
  int  szreal; /* sizeof real */
  int  pad0;   /* padding */
  real phi; /* cis is zero, clockwise positive */
  real cos; /* cos(m, n) */
  real sgn; /* (0, pi) is 1.0, otherwise -1.0 */

  real g2;
  real g[4][3]; /* gradient for each particle */

  real div; /* the divergence */
  real d4ij, d4ik, d4jj, d4jk, d4jl, d4kk, d4kl;

  unsigned int flags; /* a copy of flags used */
  int t1, t2, t3; /* GROMACS shift indices */
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

/* compute the dihedral angle, gradient g and divergence
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
      } /* do divergence */

    } else { /* clear the gradients */
      int j;
      for (j = 0; j < 4; j++)
        rv3_zero(dih->g[j]);
    }
  }
  return phi;
}

#endif

