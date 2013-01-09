/* perturbation approach to measure pressure */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_RPT
#include "zcom.h"
#include "vmove.h"

int n = 108;
int d = 3;
real rho = 0.6f; /* initial density, help set the rdf boundary */
real rcdef = 1000.0f; /* since we change volume, use half-box cutoff */
real tp = 2.5f;
real pressure = 1.0f;
real escl = 1.0f;
real amp = 0.5f;
int nsteps = 2000000;
int usesq = 1;
int nstvmove = 10; /* number of steps for volume moves */
real ra = .95f, rb = 1.5f;
av_t avU, avp, avrho;

int main(void)
{
  int t, acc = 0, vacc = 0, vtot = 0, isrun = 0;
  lj_t *lj;
  real u, p, rho1, bp, pr, pr1 = 0, bet;
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

  bet = 1.0f/tp;
  if (usesq) bet *= escl;

  for (t = 1; t <= nsteps; t++) {
    acc += lj_metro3d(lj, amp, bet);
  
    isrun = (t > nsteps / 2);

    if (t % 10 == 0) {
      if (usesq) {
        /* fixed move size */
/*
        sq_vmove(lj, -0.01, tp/escl, rpt, 1);
*/
        sq_vmove(lj, .01, tp/escl, rpt, 0);
        //sq_vmove0(lj, -0.01, tp/escl, rpt, 1);
      } else {
        lj_vmove(lj, 0.01, tp, rpt, 1);
      }
    }

    if (nstvmove && t % nstvmove == 0) {
      vtot += 1;
      if (usesq) {
        vacc += lj_mcpsq(lj, 0.02, tp/escl, pressure/escl, 0, 1e300, 0, 0);
      } else {
        vacc += lj_mcplj(lj, 0.07, tp, pressure, 0, 1e300, 0, 0);
      }
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
  bp = rpt_bets(rpt, 1);
  pr1 = bp * tp;
  printf("a %g, b %g, escl %g, erg %g, p %g (bp %g, pp %g, %g), rho %g, "
      "acc %.2f%%, vacc %.2f%%, dof/d %g\n",
      ra, rb, escl,
      u, p, bp, pr, pr1, rho1, 100.*acc/nsteps, 100.*vacc/(1e-6 + vtot),
      1.*lj->dof/lj->d);
  hs_save(hsvol, "volmc.his", HIST_NOZEROES);
  hs_close(hsvol);
  lj_close(lj);
  rpt_wdist(rpt, "dv.his");
  rpt_close(rpt);
  return 0;
}

