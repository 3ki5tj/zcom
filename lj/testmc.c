#include "lj.c"
#include "include/av.h"

int n = 108;
int d = 3;
real rho = 0.6f; /* initial density, help set the rdf boundary */
real rcdef = 1000.0f; /* since we change volume, use half-box cutoff */
real tp = 1.5f;
real pressure = 1.0f;
real amp = 0.3f;
int nsteps = 10000000, nstrdf = 1000;
int usesq = 0;
real ra = 1.0f, rb = 1.5f;
av_t avU, avp, avrho;

int main(void)
{
  int t, acc = 0, vacc = 0, vtot = 0, isrun = 0;
  lj_t *lj;
  real u, p, rho1;
  hist_t *hsvol;

  lj = lj_open(n, d, rho, rcdef);
  /* since this is MC, the dof should be n*d
   * but may we pretend it has the same dof as MD? */
  lj->dof = n * d;
  if (usesq) lj_initsq(lj, ra, rb);
  lj_rdfopen(lj, 0.01, lj->l * 2.f);

  hsvol = hs_open(1, 0, 5.*lj->n/rho, 1.f);

  for (t = 1; t <= nsteps; t++) {
    acc += lj_metro3d(lj, amp, 1.0f/tp);
   
    isrun = (t > nsteps / 2);
    
    if (t % 10 == 0) {
      vtot += 1;
      vacc += lj_mcp(lj, 0.05, tp, pressure, 0, 1e300, 0, 0);
      if (isrun) {
        av_add(&avrho, lj->n / lj->vol);
        hs_add1ez(hsvol, lj->vol, HIST_VERBOSE);
      }
    }

    if (t % 10000 == 0) { /* refresh energy regularly, just in case */
      lj_energy(lj);
    }

    if (isrun) {
      av_add(&avU, lj->epot);
      av_add(&avp, lj_calcp(lj, tp));
      if (t % nstrdf == 0) lj_rdfadd(lj);
    }
  }
  u = av_getave(&avU)/lj->n;
  p = av_getave(&avp);
  rho1 = av_getave(&avrho);
  printf("erg %g, p %g, rho %g, acc %.2f%%, vacc %.2f%%, rdfnfr %d, dof/d %g\n",
      u, p, rho1, 100.*acc/nsteps, 100.*vacc/(1e-6 + vtot), lj->rdfnfr,
      1.*lj->dof/lj->d);
  hs_save(hsvol, "volmc.his", HIST_NOZEROES);
  hs_close(hsvol);
  lj_rdfsave(lj, "rdfmc.dat", HIST_NOZEROES);
  lj_close(lj);
  return 0;
}

