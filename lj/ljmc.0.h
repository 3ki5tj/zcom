#ifndef LJMC_H__
#define LJMC_H__



/* Lennard-Jones system: Monte Carlo routines */



/* randomly displace particle i with random amplitude */
INLINE int lj_randmv3d(lj_t *lj, real *xi, real amp)
{
  int i, d;

  i = (int) (rand01() * lj->n);
  amp /= lj->l;
  rv3_copy(xi, lj->x + i*3);
  for (d = 0; d < 3; d++) /* displacement */
    xi[d] += (real) randunif(-amp, amp);
  return i;
}



/* compute energy data for a 3D pair with a square well potential */
INLINE int lj_pairsq3d(const real *xi, const real *xj, real l,
    real ra2, real rb2, real *pdr2, int inf)
{
  real dx[3], dr2;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (pdr2) *pdr2 = dr2;
  if (dr2 < ra2) return -inf;
  else if (dr2 < rb2) return 1;
  else return 0;
}



/* return the energy change (square well) from displacing x[i] to xi */
INLINE int lj_depotsq3d(lj_t *lj, int i, const real *xi, real *rm)
{
  int j, n = lj->n, npr = 0, inf = lj->esqinf, recalc = 0;
  real l = lj->l, ra2 = lj->ra2, rb2 = lj->rb2;
  real r2o, r2n, rm2o = 0, rm2 = 0;
  rv3_t *x = (rv3_t *) lj->x;
  const real tol = 1e-5;

  if (rm) rm2o = rm2 = (*rm) * (*rm);
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    npr -= lj_pairsq3d(x[i], x[j], l, ra2, rb2, &r2o, inf);
    npr += lj_pairsq3d(xi,   x[j], l, ra2, rb2, &r2n, inf);
    if (fabs(r2o - rm2o) < tol) { /* need to re-compute rmin */
      recalc |= 1;
    }
    if (r2n < rm2) { /* new rmin is found */
      recalc |= 2; /* no need to recalc */
      rm2 = r2n;
    }
  }

  /* in order to compute the minimal distance,
   * we need to occasionally recompute the entire system */
  if (recalc == 1) { /* 0, 2, 3 are safe */
    rv3_t xio;
    rv3_copy(xio, x[i]);
    rv3_copy(x[i], xi); /* apply xi */
    lj_energysq3d(lj, x, rm);
    rv3_copy(x[i], xio); /* recover */
  } else {
    if (rm) *rm = (real) sqrt(rm2);
  }

  /* hard sphere, no collision */
  if (fabs(ra2 - rb2) < 2e-6 && npr > -inf/10 && npr < inf/10) {
    npr = 0; /* number of neighbors */
  }
  return -npr; /* increased number of pairs == decreased energy */
}



/* commit a particle displacement for a square well potential */
INLINE void lj_commitsq3d(lj_t *lj, int i, const real *xi, int du)
{
  rv3_copy(lj->x + i*3, xi);
  lj->iepot += du;
  lj->epot += du;
}



/* Metropolis for a square well */
INLINE int lj_metrosq3d(lj_t *lj, real amp, real bet)
{
  int i, du;
  real xi[3], rm;

  i = lj_randmv3d(lj, xi, amp);
  rm = lj->rmin;
  du = lj_depotsq3d(lj, i, xi, &rm);
  /* patch for bet = 0 */
  if (bet >= 0 && du > lj->esqinf/2) return 0;
  if (metroacc1(du, bet)) {
    lj_commitsq3d(lj, i, xi, du);
    lj->rmin = rm;
    return 1;
  }
  return 0;
}



/* compute energy data for a particle pair, with switched potential  */
INLINE int lj_pairsw3d(lj_t *lj, real *xi, real *xj, real *u, real *vir)
{
  real dx[3], dr2, dr, fscal, psi, ksi;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, lj->l);
  if (dr2 < lj->rc2) {
    dr = (real) sqrt(dr2);
    *u = lj_potsw(lj, dr, &fscal, &psi, &ksi);
    *vir = fscal * dr2; /* f.r */
    return 1;
  } else return 0;
}



/* return the energy change from displacing x[i] to xi */
INLINE real lj_depotsw3d(lj_t *lj, int i, real *xi, real *vir)
{
  int j, n = lj->n;
  real u = 0.f, du = 0.f, dvir = 0.f;
  rv3_t *x = (rv3_t *) lj->x;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    if (lj_pairsw3d(lj, x[i], x[j], &du, &dvir)) {
      u -= du;
      *vir -= dvir;
    }
    if (lj_pairsw3d(lj, xi, x[j], &du, &dvir)) {
      u += du;
      *vir += dvir;
    }
  }
  return u;
}



/* commit a particle displacement
 * like energysw3d, it does not set pair data, lj->pr
 * call lj_forcesw3d() if it is needed */
INLINE void lj_commitsw3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  rv3_copy(lj->x + i*3, xi);
  lj->epot0 += du;
  lj->epot += du;
  lj->vir += dvir;
}



/* Metropolis algorithm */
INLINE int lj_metrosw3d(lj_t *lj, real amp, real bet)
{
  int i;
  real xi[3], du, dvir;

  i = lj_randmv3d(lj, xi, amp);
  du = lj_depotsw3d(lj, i, xi, &dvir);
  if (metroacc1(du, bet)) {
    lj_commitsw3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}



/* compute energy data for a 3D Lennard-Jones pair */
INLINE int lj_pairlj3d(real *xi, real *xj, real l, real rc2,
    real *u, real *vir)
{
  real dx[3], dr2, invr2, invr6;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (dr2 < rc2) {
    invr2 = 1.0f / dr2;
    invr6 = invr2 * invr2 * invr2;
    *vir = invr6 * (48.f * invr6 - 24.f); /* f.r */
    *u  = 4.f * invr6 * (invr6 - 1.f);
    return 1;
  } else return 0;
}



/* return the energy change from displacing x[i] to xi */
INLINE real lj_depotlj3d(lj_t *lj, int i, real *xi, real *vir)
{
  int j, n = lj->n;
  real l = lj->l, rc2 = lj->rc2, u = 0.f, du = 0.f, dvir = 0.f;
  rv3_t *x = (rv3_t *) lj->x;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    if (lj_pairlj3d(x[i], x[j], l, rc2, &du, &dvir)) {
      u -= du;
      *vir -= dvir;
    }
    if (lj_pairlj3d(xi, x[j], l, rc2, &du, &dvir)) {
      u += du;
      *vir += dvir;
    }
  }
  return u;
}



/* commit a particle displacement */
INLINE void lj_commitlj3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  rv3_copy(lj->x + i*3, xi);
  lj->epot0 += du;
  lj->epot += du;
  lj->vir += dvir;
}



INLINE int lj_metrolj3d(lj_t *lj, real amp, real bet)
{
  int i;
  real xi[3], du = 0.f, dvir = 0.f;

  i = lj_randmv3d(lj, xi, amp);
  du = lj_depotlj3d(lj, i, xi, &dvir);
  if (metroacc1(du, bet)) {
    lj_commitlj3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}



/* return the pair energy between two particles at xi and xj */
INLINE int lj_pair3d(lj_t *lj, real *xi, real *xj, real *u, real *vir)
{
  if (lj->usesq) return lj_pairsq3d(xi, xj, lj->l, lj->ra2, lj->rb2, NULL, lj->esqinf);
  if (lj->usesw) return lj_pairsw3d(lj, xi, xj, u, vir);
  return lj_pairlj3d(xi, xj, lj->l, lj->rc2, u, vir);
}



/* return the pair energy between two particles at xi and xj */
INLINE int lj_pair(lj_t *lj, real *xi, real *xj, real *u, real *vir)
{
  return lj->d == 2 ?  lj_pair2d(lj, xi, xj, u, vir) : lj_pair3d(lj, xi, xj, u, vir);
}



/* return the energy change from displacing x[i] to xi */
INLINE real lj_depot3d(lj_t *lj, int i, real *xi, real *vir, real *rmin)
{
  if (lj->usesq) return (real) lj_depotsq3d(lj, i, xi, rmin);
  if (lj->usesw) return lj_depotsw3d(lj, i, xi, vir);
  return lj_depotlj3d(lj, i, xi, vir);
}



/* return the energy change from displacing x[i] to xi */
INLINE real lj_depot(lj_t *lj, int i, real *xi, real *vir, real *rmin)
{
  return lj->d == 2 ?  lj_depot2d(lj, i, xi, vir, rmin)
    : lj_depot3d(lj, i, xi, vir, rmin);
}



/* this is defined as a macro for du is `int' in sq case, but real in other cases */
#define lj_commit3d(lj, i, xi, du, dvir) { \
  (lj->usesq) ? lj_commitsq3d(lj, i, xi, du) : \
  (lj->usesw) ? lj_commitsw3d(lj, i, xi, du, dvir) : \
                lj_commitlj3d(lj, i, xi, du, dvir); }

/* commit a move */
#define  lj_commit(lj, i, xi, du, dvir) \
  (lj->d == 2 ? lj_commit2d(lj, i, xi, du, dvir) \
              : lj_commit3d(lj, i, xi, du, dvir); )

INLINE int lj_metro3d(lj_t *lj, real amp, real bet)
{
  if (lj->usesq) return lj_metrosq3d(lj, amp, bet);
  if (lj->usesw) return lj_metrosw3d(lj, amp, bet);
  return lj_metrolj3d(lj, amp, bet);
}



/* Metropolis algorithm */
INLINE int lj_metro(lj_t *lj, real amp, real bet)
{ return lj->d == 2 ? lj_metro2d(lj, amp, bet) : lj_metro3d(lj, amp, bet); }



/* return the energy change of locally displacing a single atom */
INLINE real lj_dupertl3d(lj_t *lj, real amp)
{
  real dvir, xi[3], rmin;
  int i;

  i = lj_randmv3d(lj, xi, amp);
  return lj_depot3d(lj, i, xi, &dvir, &rmin);
}



INLINE real lj_dupertl(lj_t *lj, real amp)
{ return lj->d == 2 ? lj_dupertl2d(lj, amp) : lj_dupertl3d(lj, amp); }



/* return the energy change by random displacements of all atoms */
INLINE real lj_dupertg3d(lj_t *lj, real amp)
{
  int i, d, iep;
  rv3_t *nx;
  real du, vir, rmin, ep0, eps, lap;

  xnew(nx, lj->n);
  amp /= lj->l; /* convert to the reduced unit */
  for (i = 0; i < lj->n; i++)
    for (d = 0; d < 3; d++)
      nx[i][d] = lj->x[i*3 + d] + (real) randunif(-amp, amp);
  du = lj_energyx3d(lj, nx, &vir, &iep, &rmin, &ep0, &eps, &lap) - lj->epot;
  free(nx);
  return du;
}



INLINE real lj_dupertg(lj_t *lj, real amp)
{ return lj->d == 2 ? lj_dupertg2d(lj, amp) : lj_dupertg3d(lj, amp); }



/* return the energy caused by inserting a random atom
   the tail correction is not applied */
INLINE real lj_duinsert3d(lj_t *lj, real *xt)
{
  int j, n = lj->n;
  real xt0[3], u, du, dvir;

  if (xt == NULL)
    for (xt = xt0, j = 0; j < 3; j++)
      xt[j] = (real) rand01();
  for (u = 0.f, j = 0; j < n; j++) /* pair energy */
    if (lj_pair(lj, xt, lj->x + 3*j, &du, &dvir))
      u += du;
  return u;
}



INLINE real lj_duinsert(lj_t *lj, real *xt)
{ return lj->d == 2 ? lj_duinsert2d(lj, xt) : lj_duinsert3d(lj, xt); }

#endif
