/* molecular dynamics in the isothermal-isobaric ensemble */
#include "lj.h"
#include "ljmd.h"
#include "ljeos.h"
#include "include/av.h"

const char *fnpos = "lj.pos";
int initload = 1;

int N = 108;
real rho = 0.7f;
real rc = 1000.0f; /* use half-box cut-off */
real tp = 1.5f;
real pressure = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;
int nequil = 20000;
int nsteps = 100000;
int barostat = 4;

int usesw = 0;
real rs = 2.0f;



int main(void)
{
  lj_t *lj;
  int t, vtot = 0, vacc = 0, isrun;
  real u, k, p, rhoav, bc = 0.f;
  static av_t avU, avK, avp, avbc, avrho;
  hist_t *hsvol;

  lj = lj_open(N, 3, rho, rc);
  if (usesw) {
    lj_initsw(lj, rs);
    printf("rc %g, rs %g, box %g\n", rc, rs, lj->l*.5f);
  }
  if (initload && 0 != lj_readpos(lj, lj->x, lj->v, fnpos, LJ_LOADBOX)) {
    fprintf(stderr, "error loading previous coordinates from %s\n", fnpos);
  }

  hsvol = hist_open(1, 0, 5.*lj->n/lj->rho, 2.f);

  for (t = 1; t <= nequil + nsteps; t++) {
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    lj_vrescalex(lj, tp, thermdt);

    isrun = (t > nequil);
    if ( t % 5 == 0 ) {
      /* different barostats */
      if ( barostat == 1 ) {
        vacc += lj_mctp(lj, 0.05, tp, pressure, 0, 1e300, 0, 0);
      } else if ( barostat == 2 ) {
        vacc += lj_mcp(lj, 0.05, tp, pressure, 0, 1e300, 0, 0);
      } else if ( barostat == 3 ) {
        lj_langtp0(lj, 1e-5, tp, pressure, 0);
      } else if ( barostat == 4 ) {
        lj_langp0(lj, 1e-5, tp, pressure, 0);
      }
      vtot++;
      if (isrun) {
        av_add(&avrho, lj->n / lj->vol);
        hist_add1ez(hsvol, lj->vol, HIST_VERBOSE);
      }
    }

    if (isrun) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj_calcp(lj, tp)); /* or lj_calcpk(lj) */
      if (usesw) {
        bc = lj_bconfsw3d(lj, NULL);
        av_add(&avbc, bc);
      }
    }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  p = av_getave(&avp);
  bc = av_getave(&avbc);
  rhoav = av_getave(&avrho);
  printf("T %g, U/N %6.3f, K/N %6.3f, p %6.3f, bc %6.3f, rho %6.3f, vacc %g%%\n",
      tp, u, k, p, bc, rhoav, 100. * vacc / vtot);
  u = ljeos3d_getx(rhoav, tp, &p, NULL, NULL, LJEOS_PVEhBHKN);
  printf("ref., u %6.3f, p %6.3f\n", u, p);
  hist_save(hsvol, "volmd.his", HIST_NOZEROES);
  hist_close(hsvol);
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  return 0;
}
