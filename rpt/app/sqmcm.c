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
real rc = 2.5f; /* cutoff distance for the Lennard-Jones potential */
real tp = 1.0f;
real amp = 0.2f;
int nequil = 10000;
int nsteps = 100000;
int nevery = 10;  /* compute temperatures every this number of steps */
int calcep = 0;  /* whether to compute LJ energy or not */
real ampp = 0.1f;
int dumax = 1000; /* histogram dimension */
real sqescl = 1.0;
char *fnehis = "ehsqmc.dat"; /* energy-increment distribution */
char *fnehislj = "ehljmc.dat"; /* from the Lennard-Jones potential */
int baselj = 0; /* run the simulation on the Lennard-Jones system, doesn't work */

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &N,      "number of particles");
  argopt_add(ao, "-T", "%r", &tp,     "temperature");
  argopt_add(ao, "-r", "%r", &rho,    "density");
  argopt_add(ao, "-a", "%r", &ra,     "closest disance of the square-well potential");
  argopt_add(ao, "-b", "%r", &rb,     "the cutoff distance of the square-well potential");  
  argopt_add(ao, "-c", "%r", &rc,     "the cutoff distance of the Lennard-Jones potential");  
  argopt_add(ao, "-0", "%d", &nequil, "number of equilibration");
  argopt_add(ao, "-1", "%d", &nsteps, "number of simulation steps");
  argopt_add(ao, "-e", "%d", &nevery, "interval of computing temperatures");  
  argopt_add(ao, "-m", "%r", &amp,    "amplitude of a MC move");
  argopt_add(ao, "-M", "%r", &ampp,   "amplitude of a perturbation");
  argopt_add(ao, "-U", "%d", &dumax,  "maximal energy change for histogram");
  argopt_add(ao, "-V", "%b", &calcep, "whether to computer Lennard-Jones energy and pressure");
  argopt_add(ao, "-S", "%r", &sqescl, "energy unit of the square well potential");
  argopt_add(ao, "-I", "%b", &baselj, "run simulation with the Lennard-Jones potential");
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehislj,"output file for that using the Lennard-Jones potential");
/*
  argopt_add(ao, "-g", NULL, &fnlog, "log file");
*/
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  die_if (rb <= ra,
    "rb %g is less than ra %g, simulating hard balls\n", rb, ra);
  argopt_dump(ao);
  argopt_close(ao);
}

/* do Monte Carlo simulation and compute the perturbation temperature */
static void domc(lj_t *lj)
{
  int t, acc = 0;
  int id, idu;
  double epslj, vir, bp0, bp1, bpi, bplj0, bplj1, bplji;
  double Ulj, Uljref, plj, pljref, bet;
  double sc = 1e-30, sdU = 0, sdlnz = -1e10, Usq = 0, dU = 0, dS = 0;
  static av_t avUsq[1], aveps[1], avUlj[1], avplj[1], avepslj[1];
  rpti_t *rpt;
  rpt_t *rptlj;

  rpt = rpti_open(-dumax, dumax, 1, RPTI_HASINF);
  rptlj = rpt_open(-1000.0, 1000.0, 0.001);
  
  bet = (lj->usesq ? (sqescl/tp) : 1.0/tp);
  for (t = 0; t < nequil; t++) { /* warm up */
    lj_metro3d(lj, amp, bet);
  }

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    acc += lj_metro3d(lj, amp, bet);

    if (t % nevery == 0) {
      real xi[3];

      id = lj_randmv3d(lj, xi, ampp);
      /* compute the increment of the square-well energy */
      idu = lj_depotsq3d(lj, id, xi);
      if (idu > -1000) { /* idu may be negative infinity for the inverse mode */
        rpti_add(rpt, idu);
      }
     
      /* compute the increment of the LJ energy */ 
      epslj = lj_depotlj3d(lj, id, xi, &vir);
      rpt_add(rptlj, epslj);
      
      if (calcep) {
        if (!lj->usesq) {
          Ulj = lj->epot;
          Usq = sqescl * lj_energysq3d(lj, (rv3_t *) lj->x);
          av_add(avUsq, Usq);
          if (Usq < 1e4) {
            sdU += dU = (Usq - Ulj)/tp;
            //printf("Ulj %g, Usq %g, dU %g\n", Ulj, Usq, dU); getchar();
            sdlnz = lnadd(sdlnz, -dU); /* ln(Zsq/Zlj) */
          }
        } else {
          Usq = sqescl * lj->epot;
          Ulj = lj_energylj3d(lj, (rv3_t *) lj->x, &lj->vir, NULL, NULL);
          av_add(avUlj, Ulj);
          av_add(avplj, lj_calcp(lj, tp));
          sdU += dU = (Ulj - Usq)/tp;
          //printf("Ulj %g, Usq %g, dlnz %g\n", Ulj, Usq, sdlnz);
          sdlnz = lnadd(sdlnz, -dU); /* ln(Zlj/Zsq) */
        }
        sc += 1;
      }
      if (idu < 100) {
        av_add(aveps, sqescl * idu);
        av_add(avepslj, epslj);
      }
    }
    if (!lj->usesq) {
      av_add(avUlj, lj->epot);
    } else {
      av_add(avUsq, lj->epot * sqescl);
    }
    //if (acc) {   printf("%d, %d, %g", idu, lj->iepot, lj->epot); getchar(); }
  }
  Usq = av_getave(avUsq)/lj->n;
  Ulj = av_getave(avUlj)/lj->n;
  plj = av_getave(avplj);
  Uljref = lj_eos3d(rho, tp, &pljref, 0, 0);
  bp1 = rpti_bet1(rpt, &bp0)/sqescl;
  bp0 /= sqescl;
  bpi = rpti_bet(rpt, 0)/sqescl;
  bplj0 = rpt_bet0(rptlj);
  bplj1 = rpt_bet1(rptlj);
  bplji = rpt_bet(rptlj, 0);

  /* compute the relative entropy */
  if (sc > 0.1) {
    dU = sdU / sc;
    sdlnz -= log(sc);
    dS = dU + sdlnz;
  }
  printf("epot %g, acc %.3f, desq %5.3f(%5.3f), delj %5.3f(%5.3f), "
         "bp0 %.6f, bp1 %.6f, bpi %.6f, "
         "Ulj %.6f, Uljref %.6f, plj %.6f, pljref %.6f, "
         "bplj0 %.6f, bplj1 %.6f, bplji %.6f, "
         "dSvar %g, dU %g dlnZ %g\n",
         Usq, 1.*acc/nsteps, av_getave(aveps), av_getdev(aveps),
         av_getave(avepslj), av_getdev(avepslj),
         bp0, bp1, bpi,
         Ulj, Uljref, plj, pljref, bplj0, bplj1, bplji,
         dS, dU, sdlnz);
     
  rpti_wdist(rpt, fnehis);
  rpt_wdist(rptlj, fnehislj);
  rpti_close(rpt);
  rpt_close(rptlj);     
}

int main(int argc, char **argv)
{
  lj_t *lj;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rc);
  lj_initsq(lj, ra, rb);
  if (baselj) {
    lj->usesq = 0;
    lj_energy(lj);
  }

  domc(lj);

  lj_close(lj);
  return 0;
}
