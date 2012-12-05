/* random perturbation temperature with Monte Carlo simulation
   applied to square-well potential,
   try to match the Lennard-Jones potential */
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
real amp = 0.2f;
int nequil = 10000;
int nsteps = 100000;
int nevery = 10;  /* compute temperatures every this number of steps */
int calcUlj = 0;  /* whether to compute LJ energy or not */
real scaleU = 1.0; /* scaling for hard-sphere potential */
real ampp = 0.1f;
int dumax = 1000; /* histogram dimension */
char *fnehis = "ehsqmc.dat"; /* energy-increment distribution */
char *fnehislj = "ehljmc.dat"; /* from the Lennard-Jones potential */

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
  argopt_add(ao, "-V", "%b", &calcUlj, "whether to computer Lennard-Jones energy and pressure");
  argopt_add(ao, "-S", "%r", &scaleU,  "potential scaling for hard-sphere potential");
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehislj,"output file for that using the Lennard-Jones potential");
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
static void domc(lj_t *lj)
{
  int t, acc = 0;
  int id, idu;
  double u, eps, epslj, vir, bp0, bp1, bpi, bplj0, bplj1, bplji;
  double Ulj, Uljref, plj, pljref;
  static av_t avU[1], aveps[1], avUlj[1], avplj[1], avepslj[1];
  rpti_t *rpt;
  rpt_t *rptlj;

  rpt = rpti_open(-dumax, dumax, 1, RPTI_HASINF);
  rptlj = rpt_open(-100.0, 100.0, 0.001);
  
  for (t = 0; t < nequil; t++) /* warm up */
    lj_metrosq3d(lj, amp, scaleU/tp);

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    acc += lj_metrosq3d(lj, amp, scaleU/tp);

    if (t % nevery == 0) {
      real xi[3];
      id = lj_randmv3d(lj, xi, ampp);
      idu = lj_depotsq3d(lj, id, xi);
      rpti_add(rpt, idu);
      
      epslj = lj_depotlj3d(lj, id, xi, &vir);
      rpt_add(rptlj, epslj);
      
      //if (!metroacc1(idu, 1.0f/tp)) idu = 0;
      //rpti_add(rptd, idu);
      
      if (calcUlj) {
        Ulj = lj_energylj3d(lj, (rv3_t *) lj->x, &lj->vir, NULL, NULL);
        av_add(avUlj, Ulj);
        av_add(avplj, lj_calcp(lj, tp));
      }
      av_add(aveps, (double) idu);
      av_add(avepslj, epslj);
    }
    av_add(avU, lj->epot);
    //if (acc) {   printf("%d, %d, %g", idu, lj->iepot, lj->epot); getchar(); }
  }
  u = av_getave(avU)/lj->n;
  Ulj = av_getave(avUlj)/lj->n;
  plj = av_getave(avplj);
  Uljref = lj_eos3d(rho, tp, &pljref, 0, 0);
  eps = av_getdev(aveps);
  bp1 = rpti_bet1(rpt, &bp0)/scaleU; bp0 /= scaleU;
  bpi = rpti_bet(rpt, 0)/scaleU;
  epslj = av_getdev(avepslj);
  bplj0 = rpt_bet0(rptlj);
  bplj1 = rpt_bet1(rptlj);
  bplji = rpt_bet(rptlj, 0);

  printf("epot %g, acc %g, epsdev %g, bp0 %.6f, bp1 %.6f, bpi %.6f, "
         "Ulj %.6f, Uljref %.6f, plj %.6f, pljref %.6f, "
         "epsljdev %g, bplj0 %.6f, bplj1 %.6f, bplji %.6f\n",
     u, 1.*acc/nsteps, eps, bp0, bp1, bpi,
     Ulj, Uljref, plj, pljref,
     epslj, bplj0, bplj1, bplji);
     
  rpti_wdist(rpt, fnehis);
  rpt_wdist(rptlj, fnehislj);
  rpti_close(rpt);
  rpt_close(rptlj);     
}

int main(int argc, char **argv)
{
  lj_t *lj;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rb);
  lj_initsq(lj, ra, rb);

  domc(lj);

  lj_close(lj);
  return 0;
}
