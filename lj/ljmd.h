#ifndef LJMD_H__
#define LJMD_H__
/* Lennard-Jones system: molecular dynamics routines */


/* velocity scaling for regular (no thermostat) MD during equilibration
 * `tp' is the target temperature
 * `ekt' is the observed average kinetic energy over several steps */
#define lj_vscale(lj, tp, ekt) \
  md_vscale(lj->v, lj->n * lj->d, lj->dof, tp, ekt, &lj->ekin, &lj->tkin)

#define lj_vrescale(lj, tp, thermdt) \
  md_vrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin)

#define lj_vrescalex(lj, tp, thermdt) \
  md_vrescalex(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin);

#define lj_mcvrescale(lj, tp, thermdt) \
  md_mcvrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin);



#define lj_vv(lj, dt) lj_vvx(lj, 1.f, dt)

/* velocity Verlet */
INLINE void lj_vvx(lj_t *lj, real fscal, real dt)
{
  int i, nd = lj->n*lj->d;
  real dtl = dt/lj->l, dthf = dt * .5f * fscal;

  for (i = 0; i < nd; i++) { /* VV part 1 */
    lj->v[i] += lj->f[i] * dthf;
    lj->x[i] += lj->v[i] * dtl;
  }
  lj_force(lj); /* calculate the new force */
  for (i = 0; i < nd; i++) /* VV part 2 */
    lj->v[i] += lj->f[i] * dthf;

  lj->ekin = md_ekin(lj->v, nd, lj->dof, &lj->tkin);
  lj->t += dt;
}




/* Nose-Hoover thermostat/barostat
 * set cutoff to half of the box */
#define lj_hoovertp(lj, dt, tp, pext, zeta, eta, Q, W, ensx) \
  md_hoovertp(lj->v, lj->n, lj->d, lj->dof, dt, tp, pext, zeta, eta, \
      Q, W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin)



/* Nose-Hoover chain thermostat/barostat
 * set cutoff to half of the box */
#define lj_nhchaintp(lj, dt, tp, pext, zeta, eta, Q, M, W, ensx) \
  md_nhchaintp(lj->v, lj->n, lj->d, lj->dof, dt, tp, pext, zeta, eta, \
      Q, M, W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin)



/* Langevin barostat, with kinetic-energy scaling */
#define lj_langtp(lj, dt, tp, pext, zeta, eta, W, ensx) \
  md_langtp(lj->v, lj->n, lj->d, dt, tp, pext, zeta, eta, \
      W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin)



/* position Langevin barostat, with kinetic-energy scaling */
INLINE void lj_langtp0(lj_t *lj, real barodt, real tp, real pext, int ensx)
{
  md_langtp0(lj->v, lj->n, lj->d, barodt, tp, pext, &lj->vol,
      lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin);
  lj_setrho(lj, lj->n/lj->vol);
  lj_force(lj);
}



/* old interface */
#define lj_lgvvolmove(lj, barodt, tp, p) \
  lj_langp0(lj, barodt, tp, p, 0)

/* Langevin barostat, with coordinates only, barodt ~ 1e-5 for n = 108 */
INLINE void lj_langp0(lj_t *lj, real barodt, real tp, real pext, int ensx)
{
  md_langp0(lj->dof, lj->d, barodt, tp, pext, &lj->vol, lj->vir, lj->p_tail, ensx);
  lj_setrho(lj, lj->n/lj->vol);
  lj_force(lj);
}

/* velocity Verlet with the scaling step in the Nose-Hoover barostat */
INLINE void lj_vv_hoovertp(lj_t *lj, real dt, real eta)
{
  int i, nd = lj->n*lj->d;
  real dt2 = dt * .5f, xp;

  for (i = 0; i < nd; i++) /* VV part 1 */
    lj->v[i] += lj->f[i] * dt2;

  /* position update with scaling */
  md_hoovertpdr(lj->x, lj->v, nd, &xp, lj->l, eta, dt);
  lj->l *= xp;
  lj_setrho(lj, lj->rho/(xp*xp*xp));
  lj_force(lj); /* calculate the new force */

  for (i = 0; i < nd; i++) /* VV part 2 */
    lj->v[i] += lj->f[i] * dt2;

  lj->ekin = md_ekin(lj->v, nd, lj->dof, &lj->tkin);
  lj->t += dt;
}



/* Berendsen barostat: as a backup for constant pressure simulation */
INLINE void lj_pberendsen(lj_t *lj, real barodt, real tp, real pext)
{
  int i;
  real pint, vn, lo = lj->l, s, dlnv;

  pint = (lj->vir + 2.f * lj->ekin)/ (lj->d * lj->vol) + lj->p_tail;

  /* proposed change of log V */
  dlnv = (pint - pext)*lj->vol/tp*barodt;
  if (dlnv < -0.1) dlnv = -0.1; else if (dlnv > 0.1) dlnv = 0.1;
  vn = log(lj->vol) + dlnv;
  vn = exp(vn);
  lj_setrho(lj, lj->n/vn);
  s = lo/lj->l;
  for (i = 0; i < lj->d * lj->n; i++) lj->v[i] *= s;
  lj->ekin *= s*s;
  lj->tkin *= s*s;
}



/* In Monte Carlo barostats, we compute the energy directly */
#define LJ_FIXEDRC 0x4000

#define lj_mcprescale(lj, lnvamp, tp, pext, vmin, vmax, ensx) \
  lj_mctp(lj, lnvamp, tp, pext, vmin, vmax, ensx, 0)

/* Monte Carlo barostat, with kinetic-energy scaling
 * the ideal gas part computed as \sum p^2/m / V
 * the scaling is r = r*s, p = p/s;
 * set cutoff to half of the box */
INLINE int lj_mctp(lj_t *lj, real lnvamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, s, epo, bet = 1.f/tp;
  double dex;
  lj_t *lj1;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;
  if (vn < vmin || vn >= vmax)
    return 0;
  if ((flags & LJ_FIXEDRC) && ln < lj->rc * 2)
    return 0; /* box too small */

  epo = lj->epot;
  lj1 = lj_clone(lj, LJ_CPF); /* make a copy */
  lj_setrho(lj, lj->n/vn);
  lj_force(lj); /* we change force here */
  dex = bet * (lj->epot - epo + pext * (vn - vo))
      + bet * (pow(vo/vn, 2.0/d) - 1)*lj->ekin
      + d * (lnlo - lnln) * (1 - ensx);
  if (metroacc1(dex, 1.f)) { /* scale the velocities */
    s = lo/lj->l;
    for (i = 0; i < d * lj->n; i++) lj->v[i] *= s;
    lj->ekin *= s*s;
    lj->tkin *= s*s;
    acc = 1;
  } else {
    lj_copy(lj, lj1, LJ_CPF); /* restore force etc. */
  }
  lj_close(lj1);
  return acc;
}



/* Monte Carlo barostat for the square-well potential, for coordinates only
 * suppose lj->rmin has been correctly set
 * use lnvamp 0.03 ~ 0.06 for 256 system */
INLINE int lj_mcpsq(lj_t *lj, real lnvamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags)
{
  int acc = 0, i, d = lj->d, iep;
  real lnlo, lnln, vo, vn, lo, ln, rmn = 0, epo, bet = 1.f/tp;
  double dex;

  (void) flags;
  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;
  if (vn < vmin || vn >= vmax) return 0;

  /* check if there is a clash */
  rmn = lj->rmin * ln / lo;
  if (ln < lo) {
    if (ln < lj->rb * 2) return 0; /* box too small */
    if (rmn < lj->ra) return 0;
  }

  /* compute the change of the square-well energy */
  epo = lj->epot;
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  if (fabs(lj->ra - lj->rb) < 1e-6) { /* skip the energy calculation */
    iep = 0;
  } else {
    if (d == 3) {
      iep = lj_energysq3d(lj, (rv3_t *) lj->x, &rmn);
    } else {
      iep = lj_energysq2d(lj, (rv2_t *) lj->x, &rmn);
    }
  }
  dex = bet * ((real) iep - epo + pext * (vn - vo))
      - (lj->dof + (1 - ensx) * d) * (lnln - lnlo);
  if (rmn > lj->ra && metroacc1(dex, 1.0)) {
    lj->iepot = iep;
    lj->epot = iep;
    lj->rmin = rmn;
    acc = 1;
  } else {
    lj_setrho(lj, lj->n/vo);
  }
  return acc;
}



/* Monte Carlo barostat, coordinates only
 * use lnvamp 0.03 ~ 0.06 for 256 system */
INLINE int lj_mcplj(lj_t *lj, real lnvamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, epo, bet = 1.f/tp;
  double dex;
  lj_t *lj1;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;
  if (vn < vmin || vn >= vmax)
    return 0;
  if ((flags & LJ_FIXEDRC) && ln < lj->rc * 2)
    return 0; /* box too small */

  epo = lj->epot;
  lj1 = lj_clone(lj, LJ_CPF); /* save a copy */
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  lj_force(lj);
  dex = bet * (lj->epot - epo + pext * (vn - vo))
      - (lj->dof + (1 - ensx) * d) * (lnln - lnlo);
  if (metroacc1(dex, 1.0)) {
    acc = 1;
  } else {
    lj_copy(lj, lj1, LJ_CPF);
  }
  lj_close(lj1);
  return acc;
}



/* old interface */
#define lj_volmove(lj, lnlamp, tp, p) \
  lj_mcp(lj, lnlamp*lj->d, tp, p, 0, 1e300, 0, LJ_FIXEDRC)

/* Monte Carlo barostat, coordinates only */
INLINE int lj_mcp(lj_t *lj, real lnvamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags)
{
  if (lj->usesq) { /* use the specialized square-well version */
    return lj_mcpsq(lj, lnvamp, tp, pext, vmin, vmax, ensx, flags);
  } else { /* use the generic version */
    return lj_mcplj(lj, lnvamp, tp, pext, vmin, vmax, ensx, flags);
  }
}

#endif
