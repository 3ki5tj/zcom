
/* MC virtual volume moves for the square-well potential */
INLINE int sq_vmove(lj_t *lj, real lnvamp, real tp, rpt_t *rpt, int fixed)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, epo, epn, bet = 1.f/tp, rmn, rmn0 = 1e199;
  double dex, wt;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  if (fixed) {
    lnln = lnlo + lnvamp/d;
  } else {
    lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  }
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;

  /* check if there is a clash */
  if (ln < lo) {
    if (ln < lj->rb * 2) goto ERR; /* box too small */
    rmn = lj->rmin * ln / lo;
    rmn0 = rmn;
    if (rmn < lj->ra) goto ERR;
  }

  /* for a hard-sphere potential, no energy, accept it immediately */
  if (fabs(lj->ra - lj->rb) < 1e-6) {
    epn = epo = 0;
    goto ACC;
  }

  epo = lj->epot;
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  if (lj->d == 3) {
    epn = lj_energysq3d(lj, (rv3_t *) lj->x, &rmn);
  } else {
    epn = lj_energysq2d(lj, (rv2_t *) lj->x, &rmn);
  }
  if (rmn < lj->ra) { /* should never happen, user might change ra */
    fprintf(stderr, "undetected collision! rmn %g, ra %g\n", rmn, lj->ra);
    goto ERR;
  }
  lj_setrho(lj, lj->n/vo);
ACC:
  dex = bet * (epn - epo) - (lj->dof + d) * (lnln - lnlo);
  wt = exp(-dex);
  if (fabs(vn - vo) > 100) {
    fprintf(stderr, "add success vo %g, vn %g, ra %g, rmn %g, %g\n", vo, vn, lj->ra, rmn0, rmn);
    getchar();
  }  
  rpt_addw(rpt, vn - vo, wt);
  return acc;
ERR:
  if (fabs(vn - vo) > 100) {
    fprintf(stderr, "add fail vo %g, vn %g\n", vo, vn);
    getchar();
  }
  rpt_addw(rpt, vn - vo, 0);
  return 0;
}

/* back-up routine for sq_vmove */
INLINE int sq_vmove0(lj_t *lj, real lnvamp, real tp, rpt_t *rpt, int fixed)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, epo, epn, bet = 1.f/tp, rmn;
  double dex, wt;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  if (fixed) {
    lnln = lnlo + lnvamp/d;
  } else {
    lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  }
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;

  epo = lj->epot;
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  if (lj->d == 3) {
    epn = lj_energysq3d(lj, (rv3_t *) lj->x, &rmn);
  } else {
    epn = lj_energysq2d(lj, (rv2_t *) lj->x, &rmn);
  }
  dex = bet * (epn - epo) - (lj->dof + d) * (lnln - lnlo);
  wt = exp(-dex);
  rpt_addw(rpt, vn - vo, wt);
  lj_setrho(lj, lj->n/vo);
  return acc;
}


/* MC-like virtual move
 * set cutoff to half of the box */
INLINE int lj_vmove(lj_t *lj, real lnvamp, real tp, rpt_t *rpt, int fixed)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, epo, epn, bet = 1.f/tp, vir, ep0, eps;
  double dex, wt;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  if (fixed) {
    lnln = lnlo + lnvamp/d;
  } else {
    lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  }
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;

  epo = lj->epot;
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  if (lj->d == 3) {
    epn = lj_energylj3d(lj, (rv3_t *) lj->x, &vir, &ep0, &eps);
  } else {
    epn = lj_energylj2d(lj, (rv2_t *) lj->x, &vir, &ep0, &eps);
  }
  dex = bet * (epn - epo) - (lj->dof + d) * (lnln - lnlo);
  wt = exp(-dex);
  rpt_addw(rpt, vn - vo, wt);
  lj_setrho(lj, lj->n/vo);
  return acc;
}


