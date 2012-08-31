/* random perturbation temperature with Monte Carlo simulation
   applied to square-well potential 
   only use adujested move */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_RPT
#define ZCOM_ARGOPT
#include "zcom.h"

#define D 3

int N = 108;
real rho = 0.7f;
real ra = 1.0f; /* barrier distance of the square-well potential */
real rb = 1.5f; /* cutoff distance of the square-well potential */
real tp = 1.0f;
real amp = 0.04f;
int nequil = 10000;
int nsteps = 100000;
int nevery = 1;  /* compute temperatures every this number of steps */


real ampp = 0.04f;
int dumax = 1000; /* histogram dimension */
char *fnehis = "ehsqmc.dat"; /* energy-increment distribution */
char *fnehisd = NULL; /* adjusted energy-increment distribution */

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,      "number of particles");
  argopt_add(ao, "-T", "%r", &tp,     "temperature");
  argopt_add(ao, "-r", "%r", &rho,    "density");
  argopt_add(ao, "-a", "%r", &ra,     "closest disance of the square-well potential");
  argopt_add(ao, "-b", "%r", &rb,     "the cutoff distance of the square-well potential");  
  argopt_add(ao, "-0", "%d", &nequil, "number of equilibration");
  argopt_add(ao, "-1", "%d", &nsteps, "number of simulation steps");
  argopt_add(ao, "-e", "%d", &nevery, "interval of computing temperatures");  
  argopt_add(ao, "-m", "%r", &amp,    "amplitude of a MC move");
  argopt_add(ao, "-M", "%r", &ampp,   "amplitude of a perturbation");
  argopt_add(ao, "-U", "%d", &dumax,  "maximal energy change for histogram");
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehisd,"output file for the adjusted energy-increment histogram");
/*
  argopt_add(ao, "-g", NULL, &fnlog, "log file");
*/
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  if (rb <= ra) {
    printf("rb %g is less than ra %g, simulating hard balls\n", rb, ra);
    exit(1); /* do not allow it */
  }
  argopt_dump(ao);
  argopt_close(ao);
}

/* do Monte Carlo simulation and compute the perturbation temperature */
static void domc(lj_t *lj, rpti_t *rpt, rpti_t *rptd)
{
  int t, acc = 0;
  int id, idu;
  double u, du, bp0, bp1, bpi, bps0, bps1, bph0, bph1, bpd0, bpd1, bpdi;
  static av_t avU[1], avu[1];

  for (t = 0; t < nequil; t++) /* warm up */
    lj_metrosq3d(lj, amp, 1.0f/tp);

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    acc += lj_metrosq3d(lj, amp, 1.0f/tp);

    if (t % nevery == 0) {
      real xi[3];
      id = lj_randmv3d(lj, xi, ampp);
      idu = lj_depotsq3d(lj, id, xi);
      rpti_add(rpt, idu);
      
      if (!metroacc1(idu, 1.0f/tp)) idu = 0;
      rpti_add(rptd, idu);
      
      av_add(avu, (double) idu);
    }
    av_add(avU, lj->epot);
    //if (acc) {   printf("%d, %d, %g", idu, lj->iepot, lj->epot); getchar(); }
  }
  u = av_getave(avU)/lj->n;
  du = av_getdev(avu);
  bp1 = rpti_bet1(rpt, &bp0);
  bpi = rpti_bet(rpt);
  bps0 = rpti_bets(rpt, 0);
  bps1 = rpti_bets(rpt, 1);
  bph0 = rpti_beth(rpt, 0);
  bph1 = rpti_beth(rpt, 1);
  bpd1 = rpti_bet1(rptd, &bpd0);
  bpdi = rpti_bet(rptd);
  printf("epot %g, acc %g, dudev %g, bp0 %.6f, bp1 %.6f, bpi %.6f, bps0 %.6f, bps1 %.6f, "
    "bph0 %.6f, bph1 %.6f, bpd0 %.6f, bpd1 %.6f, bpdi %.6f\n",
     u, 1.*acc/nsteps, du, bp0, bp1, bpi, bps0, bps1, 
     bph0, bph1, bpd0, bpd1, bpdi);
}

int main(int argc, char **argv)
{
  lj_t *lj;
  rpti_t *rpt, *rptd;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rb);
  lj_initsq(lj, ra, rb);
  rpt = rpti_open(-dumax, dumax, 1, RPTI_HASINF);
  rptd = rpti_open(-dumax, dumax, 1, 0);
   
  domc(lj, rpt, rptd);
  
  rpti_wdist(rpt, fnehis);
  if (fnehisd) rpti_wdist(rptd, fnehisd);
  rpti_close(rpt);
  rpti_close(rptd);
  lj_close(lj);
  return 0;
}
