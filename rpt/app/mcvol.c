/* perturbation approach to measure pressure */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_RPT
#include "zcom.h"

int n = 108;
int d = 3;
real rho = 0.6f; /* initial density, help set the rdf boundary */
real rcdef = 1000.0f; /* since we change volume, use half-box cutoff */
real tp = 1.5f;
real pressure = 1.0f;
real amp = 0.3f;
int nsteps = 1000000;
int usesq = 0;
int nstvmove = 10; /* number of steps for volume moves */
real ra = 1.0f, rb = 1.5f;
av_t avU, avp, avrho;

/* MC virtual volume moves */
INLINE int sq_vmove(lj_t *lj, real lnvamp, real tp, rpt_t *rpt)
{
  int i, d = lj->d, iep;
  real lnlo, lnln, vo, vn, lo, ln, rmn = 0, epo, bet = 1.f/tp;
  double dex, wt;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
  ln = (real) exp(lnln);
  for (vn = 1, i = 0; i < d; i++) vn *= ln;

  /* check if there is a clash */
  if (ln < lo) {
    if (ln < lj->rb * 2) return 0; /* box too small */
    rmn = lj->rmin * ln / lo;
    if (rmn < lj->ra) return 0;
  }

  /* for a hard-sphere potential, no energy, accept it immediately */
  if (fabs(lj->ra - lj->rb) < 1e-6) {
    lj->rmin *= ln/lo;
    lj_setrho(lj, lj->n/vn);
    return 1;
  }

  /* compute the change of the square-well energy */
  epo = lj->epot;
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  if (d == 3) {
    iep = lj_energysq3d(lj, (rv3_t *) lj->x, &rmn);
  } else {
    iep = lj_energysq2d(lj, (rv2_t *) lj->x, &rmn);
  }

  if (rmn > lj->ra) {
    dex = bet * ((real) iep - epo) - (lj->dof + d) * (lnln - lnlo);
    wt = exp(-dex);
    rpt_addw(rpt, vn - vo, wt);
  }

  lj_setrho(lj, lj->n/vo);
  return 0;
}

/* MC-like virtual move
 * set cutoff to half of the box */
INLINE int lj_vmove(lj_t *lj, real lnvamp, real tp, rpt_t *rpt)
{
  int acc = 0, i, d = lj->d;
  real lnlo, lnln, lo, ln, vo, vn, epo, epn, bet = 1.f/tp, vir, ep0, eps;
  double dex, wt;

  vo = lj->vol;
  lo = lj->l;
  lnlo = (real) log(lo);
  lnln = (real) (lnlo + lnvamp/d * (2.f * rnd0() - 1.f));
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

int main(void)
{
  int t, acc = 0, vacc = 0, vtot = 0, isrun = 0;
  lj_t *lj;
  real u, p, rho1, bp, pr;
  hist_t *hsvol;
  rpt_t *rpt;

  lj = lj_open(n, d, rho, rcdef);
  /* since this is MC, the dof should be n*d
   * maybe we can that it were the same dof as MD? */
  lj->dof = n * d;
  if (usesq) {
    lj_initsq(lj, ra, rb);
    lj->vir = 0;
  }

  hsvol = hs_open(1, 0, 5.*lj->n/rho, 1.f);

  rpt = rpt_open(-100.0, 100.0, 0.001f);

  for (t = 1; t <= nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
   
    isrun = (t > nsteps / 2);
    
    if (t % 10 == 0) {
      if (usesq) sq_vmove(lj, 0.01, tp, rpt);
      else lj_vmove(lj, 0.01, tp, rpt);
    }

    if (nstvmove && t % nstvmove == 0) {
      vtot += 1;
      vacc += lj_mcp(lj, 0.01, tp, pressure, 0, 1e300, 0, 0);
      if (isrun) {
        av_add(&avrho, lj->n / lj->vol);
        hs_add1ez(hsvol, lj->vol, HIST_VERBOSE);
      }
    }

    if (isrun) {
      av_add(&avU, lj->epot);
      if (!lj->usesq) av_add(&avp, lj_calcp(lj, tp));
    }
  }
  u = av_getave(&avU)/lj->n;
  p = av_getave(&avp);
  rho1 = av_getave(&avrho);
  bp = rpt_betw(rpt);
  pr = bp * tp;
  printf("erg %g, p %g (pp %g), rho %g, acc %.2f%%, vacc %.2f%%, dof/d %g\n",
      u, p, pr, rho1, 100.*acc/nsteps, 100.*vacc/(1e-6 + vtot),
      1.*lj->dof/lj->d);
  hs_save(hsvol, "volmc.his", HIST_NOZEROES);
  hs_close(hsvol);
  lj_close(lj);
  return 0;
}

