/* random perturbation temperature with molecular dynamics simulation */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_RPT
#define ZCOM_ARGOPT
#include "zcom.h"

#define D 3

int N = 108;
real rho = 0.7f;
real rcdef = 2.5f; /* cut-off distance */
real rcshf = 2.0f; /* shift distance */
real tp = 1.0f;
real mddt = 0.002f; /* time step for molecular dynamics */
int mequil = 7;  /* number of times to adjust the temperature in equilibration */
int nequil = 1000;  /* number of steps of each equilibration adjustment */
int nsteps = 10000;
int nevery = 1;  /* compute temperatures every this number of steps */
int usesw = 0;  /* use the switched LJ potential */

int gpert = 0;  /* use global perturbation */
real ampp = 0.05f; /* amplitude for perturbation */
real dumax = 100.0f, dudel = 0.001f; /* histogram dimension */
char *fnehis = "eh.dat"; /* energy-increment distribution */
char *fnehisd = NULL;

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,      "number of particles");
  argopt_add(ao, "-T", "%r", &tp,     "temperature");
  argopt_add(ao, "-w", "%b", &usesw,  "use the switched LJ potential");
  argopt_add(ao, "-r", "%r", &rho,    "density");
  argopt_add(ao, "-c", "%r", &rcdef,  "cutoff distance");
  argopt_add(ao, "-s", "%r", &rcshf,  "shift distance");
  argopt_add(ao, "-j", "%d", &mequil, "number of times to adjust the kinetic energy during equilibration");
  argopt_add(ao, "-0", "%d", &nequil, "number of equilibration steps in each adjustment section");
  argopt_add(ao, "-1", "%d", &nsteps, "number of simulation steps");
  argopt_add(ao, "-e", "%d", &nevery, "interval of computing temperatures");
  argopt_add(ao, "-d", "%r", &mddt,   "time step for molecular dynamics");
  argopt_add(ao, "-M", "%r", &ampp,   "amplitude of a perturbation test");
  argopt_add(ao, "-g", "%b", &gpert,  "use global perturbation");
  argopt_add(ao, "-U", "%r", &dumax,  "maximal energy change for histogram");
  argopt_add(ao, "-u", "%r", &dudel,  "energy bin size for histogram");
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehisd, "output file for the adjusted energy-increment histogram");
/*
  argopt_add(ao, "-g", NULL, &fnlog, "log file");
*/
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  if (rcshf > rcdef) {
    rcshf = rcdef;
    if (usesw) printf("setting shift distance to cutoff %g\n", rcshf);
  }
  argopt_dump(ao);
  argopt_close(ao);
}

/* do Monte Carlo simulation and compute the perturbation temperature */
static void domc(lj_t *lj, rpt_t *rpt, rpt_t *rptd)
{
  int im, t;
  real ek, u, p, du;
  double bc, bc0, bcr, bp0, bp1, bpi, bps0, bps1, bph0, bph1, bpd0, bpd1, bpdi;
  static av_t avU[1], avp[1], avbc[1], avlap[1], avf2[1], avbc0[1], avK[1];

  for (im = 0; im < mequil; im++) {
    av_clear(avK);
    for (t = 0; t <= nequil; t++) {
      lj_vv(lj, mddt);
      lj_shiftcom(lj, lj->v);
      av_add(avK, lj->ekin);
    }
    ek = av_getave(avK);
    lj_vscale(lj, tp, ek);
    printf("t %g: v-scaling: T %g(%g), now K %g, K + U = %g\n",
        (t+1)*mddt, ek*2.f/lj->dof, avK->s, lj->ekin, lj->epot + lj->ekin);
  }

  av_clear(avK);
  for (t = 1; t <= nsteps; t++) { /* real simulation */
    lj_vv(lj, mddt);
    lj_shiftcom(lj, lj->v);
    if (t % nevery == 0) {
      lj_force3d(lj);
      bc0 = lj->lap/lj->f2;
      bc = (lj->usesw) ? lj_bconfsw3d(lj, NULL) : bc0;
      av_add(avbc0, bc0);
      av_add(avbc, bc);
      av_add(avlap, lj->lap);
      av_add(avf2, lj->f2);
      du = (gpert ? lj_dupertg3d(lj, ampp) : lj_dupertl3d(lj, ampp));
      rpt_add(rpt, du);
      //printf("t %d, bc %g, bc0 %g, du %g\n", t, bc, bc0, du);

      /* construct the adjusted change */
      if (!metroacc1(du, 1.0f/tp)) du = 0;
      rpt_add(rptd, du);
    }
    av_add(avK, lj->ekin);
    av_add(avU, lj->epot);
    av_add(avp, lj_calcp(lj, tp));
  }
  ek = av_getave(avK)/lj->n;
  u = av_getave(avU)/lj->n;
  p = av_getave(avp);
  du = av_getdev(&rpt->av);
  bc = av_getave(avbc);
  bc0 = av_getave(avbc0);
  bcr = av_getave(avlap)/av_getave(avf2);
  bp0 = rpt_bet0(rpt);
  bp1 = rpt_bet1(rpt);
  bpi = rpt_bet(rpt);
  bps0 = rpt_bets(rpt, 0);
  bps1 = rpt_bets(rpt, 1);
  bph0 = rpt_beth(rpt, 0);
  bph1 = rpt_beth(rpt, 1);
  bpd0 = 1.0f/tp + rpt_bet0(rptd);
  bpd1 = 1.0f/tp + rpt_bet1(rptd);
  bpdi = 1.0f/tp + rpt_bet(rptd);
  printf("ekin %g, epot %g, p %g, dudev %g, bc %.6f, bc0 %.6f, bcr %.6f, "
    "bp0 %.6f, bp1 %.6f, bpi %.6f, bps0 %.6f, bps1 %.6f, bph0 %.6f, bph1 %.6f,"
    "bpd0 %.6f, bpd1 %.6f, bpdi %.6f",
    ek, u, p, du, bc, bc0, bcr,
    bp0, bp1, bpi, bps0, bps1, bph0, bph1,
    bpd0, bpd1, bpdi);
}

int main(int argc, char **argv)
{
  lj_t *lj;
  rpt_t *rpt, *rptd;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rcdef);
  if (usesw) lj_initsw(lj, rcshf);
  rpt = rpt_open(-dumax, dumax, dudel);
  rptd = rpt_open(-dumax, dumax, dudel);

  domc(lj, rpt, rptd);

  rpt_wdist(rpt, fnehis);
  if (fnehisd) rpt_wdist(rptd, fnehisd);
  rpt_close(rpt);
  rpt_close(rptd);
  lj_close(lj);
  return 0;
}
