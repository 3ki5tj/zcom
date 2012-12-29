/* Bennett acceptance ratio method applied to the matching of
   the square-well potential and the Lennard-Jones potential */
#define ZCOM_PICK
#define ZCOM_LJ
#define ZCOM_RPT
#define ZCOM_ARGOPT
#include "zcom.h"

#define D 3

int N = 108;
real rho = 0.05f;
real ra = 1.0f; /* barrier distance of the square-well potential */
real rb = 1.5f; /* cutoff distance of the square-well potential */
real rc = 2.5f; /* cutoff distance for the Lennard-Jones potential */
real tp = 1.0f;
real amp = 0.8f;
int nequil = 10000;
int nsteps = 100000;
int nevery = 10;  /* compute temperatures every this number of steps */
real ampp = 0.2f;
int dumax = 1000; /* histogram dimension */
real sqescl = 1.0;
char *fnehis = "ehsqmc.dat"; /* energy-increment distribution */
char *fnehislj = "ehljmc.dat"; /* from the Lennard-Jones potential */
int verbose = 0;


int esqinf = 1000; /* infinity energy in square well potential */

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
  argopt_add(ao, "-S", "%r", &sqescl, "energy unit of the square well potential");
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehislj,"output file for that using the Lennard-Jones potential");
  argopt_add(ao, "-v", "%d", &verbose, "verbose mode");
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

/* compute dlnZ from Bennett's acceptance ratio method */ 
static double bardlnZ(hist_t *hs0, hist_t *hs1, double *dS0, double *dS1)
{
  double C = 0, Cn, *h0 = hs0->arr, *h1 = hs1->arr, x, s0, s1, sf0, sf1;
  int iter, i, i0s, i0t, i1s, i1t, n = hs0->n;

  /* search for the boundaries */
  for (i0s = 0; i0s < n; i0s++) if (h0[i0s] > 0) break;
  for (i0t = n; i0t > i0s; i0t--) if (h0[i0t - 1] > 0) break;
  for (*dS0 = 0, s0 = 0, i = i0s; i < i0t; i++) {
    if (h0[i] < 0) continue;
    x = hs0->xmin + hs0->dx * (i + .5);
    s0 += h0[i];
    *dS0 += h0[i] * x;
  }
  if (s0 > 0) *dS0 /= s0;

  for (i1s = 0; i1s < n; i1s++) if (h1[i1s] > 0) break;
  for (i1t = n; i1t > i1s; i1t--) if (h1[i1t - 1] > 0) break;
  for (*dS1 = 0, s1 = 0, i = i1s; i < i1t; i++) {
    if (h1[i] < 0) continue;
    x = hs1->xmin + hs1->dx * (i + .5);
    s1 += h1[i];
    *dS1 += h1[i] * x; 
  }
  if (s1 > 0) *dS1 /= s1;

  if (s1 <= 0) { /* no data */
    fprintf(stderr, "s0 %g, s1 %g\n", s0, s1);
    return 0;
  }
 
  /* get a rough estimate, by Z0/Z1 = < exp(U1 - U0) >_1 */
  if (s1 > 0) {
    for (sf1 = -1e100, i = i1s; i < i1t; i++) {
      if (h1[i] <= 0) continue;
      x = hs1->xmin + hs1->dx * (i + .5);
      sf1 = lnadd(sf1, -x + log(h1[i]));
    }
    sf1 -= log(s1);
    C = sf1;
  }
  
  /* get another rough estimate, by Z1/Z0 = < exp(U0 - U1) >_0 */
  if (s0 > 0) {
    for (sf0 = -1e100, i = i0s; i < i0t; i++) {
      if (h0[i] <= 0) continue;
      x = hs0->xmin + hs0->dx * (i + .5);
      sf0 = lnadd(sf0, -x + log(h0[i]));
    }
    sf0 -= log(s0);
    Cn = -sf0;
  }
  fprintf(stderr, "initial guess, s0 %g, s1 %g, dlnZ %g,%g dS0 %g, dS1 %g\n",
      s0, s1, C,Cn, *dS0, *dS1);

  if (s0 <= 0) { /* use the hs1 only */
    *dS0 -= C;
    *dS1 += C;
    return C;
  }

  for (iter = 0; iter < 1000; iter++) {
    /* compute <f(U0 - U1 + C)>_1 */
    for (sf1 = -1e100, i = i1s; i < i1t; i++) {
      if (h1[i] <= 0) continue;
      x = hs1->xmin + hs1->dx * (i + .5) + C;
      x = -lnadd(x, 0); /* 1/(1 + exp(x)) */
      sf1 = lnadd(sf1, x + log(h1[i]));
    }
    sf1 -= log(s1);

    /* compute <f(U1 - U0 - C)>_0 */
    for (sf0 = -1e100, i = i0s; i < i0t; i++) {
      if (h0[i] <= 0) continue;
      x = hs0->xmin + hs0->dx * (i + .5) - C;
      x = -lnadd(x, 0); /* 1/(1 + exp(x)) */
      sf0 = lnadd(sf0, x + log(h0[i]));
    }
    sf0 -= log(s0);
    
    Cn = C + sf1 - sf0 + log(s1/s0);
    if (verbose) {
      fprintf(stderr, "sf1 %g, sf0 %g, dlnZ %g, dlnZ' %g\n", sf1, sf0, C, Cn);
      if (verbose >= 2) getchar();
    }
    if (fabs(C - Cn) < 1e-6) break;
    C = (C + Cn)*.5; /* to prevent oscillation, but slows down convergence */
  }

  C -= log(s1/s0); /* dlnZ */
  *dS0 -= C;
  *dS1 += C;
  return C;
}

/* do Monte Carlo simulation and compute the perturbation temperature */
static void domc(lj_t *lj, lj_t *sq)
{
  int t, acclj = 0, accsq = 0;
  int id, idu;
  double epslj, vir, bp0, bp1, bpi, bplj0, bplj1, bplji;
  double dU, Ulj, Usq, UljB, UsqB, Uljref, pljB, plj, pljref;
  double dlnZ = 0, dSsq = 0, dSlj = 0;
  static av_t avUsq[1], avUlj[1], avplj[1], avUsqB[1], avUljB[1], avpljB[1], avdesq[1], avdelj[1];
  rpti_t *rptsq;
  rpt_t *rptlj;
  hist_t *hlj, *hsq; /* energy histograms */

  rptsq = rpti_open(-dumax, dumax, 1, RPTI_HASINF);
  rptlj = rpt_open(-1000.0, 1000.0, 0.001);
  hlj = hs_open(1, -20*lj->n, 20*lj->n, 0.001*lj->n);
  hsq = hs_open(1, -20*lj->n, 20*lj->n, 0.001*lj->n);
  
  for (t = 0; t < nequil; t++) { /* warm up */
    lj_metro3d(lj, amp, 1.0/tp);
    lj_metro3d(sq, amp, (sqescl/tp));  /* exp(-sqescl * de / tp) */
  }

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    acclj += lj_metro3d(lj, amp, 1.0/tp);
    accsq += lj_metro3d(sq, amp, (sqescl/tp));
    if (t % nevery == 0) {
      real xi[3];

      /* in the LJ system, compute the increment of the square-well energy */
      id = lj_randmv3d(lj, xi, ampp);
      idu = lj_depotsq3d(lj, id, xi);
      if (idu > -esqinf/2 && idu < esqinf/2) { /* idu may be negative infinity for the inverse mode */
        av_add(avdesq, idu * sqescl);
        rpti_add(rptsq, idu);
      }
      
      /* compute <Usq> in the LJ system */
      Ulj = lj->epot;
      Usq = sqescl * lj_energysq3d(lj, (rv3_t *) lj->x);
      if (Usq <= 0) /* collect configurations with no clash */
        av_add(avUsqB, Usq);
      dU = dblmin( (Usq - Ulj)/tp, hlj->xmin + hlj->dx * (hlj->n - 1) );
      hs_add1ez(hlj, dU, HIST_VERBOSE);
      
      /* in the square-well system, compute the increment of the LJ energy */ 
      id = lj_randmv3d(sq, xi, ampp);
      epslj = lj_depotlj3d(sq, id, xi, &vir);
      rpt_add(rptlj, epslj);
      av_add(avdelj, epslj);
     
      /* compute <Ulj> in the sq system */ 
      Usq = sqescl * sq->epot;
      Ulj = lj_energylj3d(sq, (rv3_t *) sq->x, &sq->vir, NULL, NULL);
      av_add(avUljB, Ulj);
      av_add(avpljB, lj_calcp(sq, tp));
      hs_add1ez(hsq, (Ulj - Usq)/tp, HIST_VERBOSE);
    }
    av_add(avUlj, lj->epot);
    av_add(avplj, lj_calcp(lj, tp));
    av_add(avUsq, sq->epot * sqescl);
  }
  Usq = av_getave(avUsq)/lj->n;
  Ulj = av_getave(avUlj)/lj->n;
  plj = av_getave(avplj);
  UsqB = av_getave(avUsqB)/lj->n;
  UljB = av_getave(avUljB)/lj->n;
  pljB = av_getave(avpljB);
  Uljref = lj_eos3d(rho, tp, &pljref, 0, 0);

  bp1 = rpti_bet1(rptsq, &bp0)/sqescl; /* we got s * bet, so divide s */
  bp0 /= sqescl;
  bpi = rpti_bet(rptsq, 0)/sqescl;
  bplj0 = rpt_bet0(rptlj);
  bplj1 = rpt_bet1(rptlj);
  bplji = rpt_bet(rptlj, 0);

  /* compute the relative entropy */
  dlnZ = bardlnZ(hlj, hsq, &dSlj, &dSsq);
  printf("accsq %.3f, acclj %.3f, desq %5.3f(%5.3f), delj %5.3f(%5.3f), "
         "Usqlj %.6f, Usq %.6f, bp0 %.6f, bp1 %.6f, bpi %.6f, "
         "Uljsq %.6f, Ulj %.6f, Uljref %.6f, pljsq %.6f, plj %.6f, pljref %.6f, "
         "bplj0 %.6f, bplj1 %.6f, bplji %.6f, "
         "dSlj %g, dSsq %g, dlnZ %g\n",
         1.*accsq/nsteps, 1.*acclj/nsteps,
         av_getave(avdesq), av_getdev(avdesq),
         av_getave(avdelj), av_getdev(avdelj),
         UsqB, Usq, bp0, bp1, bpi,
         UljB, Ulj, Uljref, pljB, plj, pljref, bplj0, bplj1, bplji,
         dSlj, dSsq, dlnZ);
     
  rpt_wdist(rptlj, fnehislj);
  rpti_wdist(rptsq, fnehis);
  rpt_close(rptlj);
  rpti_close(rptsq);
  hs_close(hlj);
  hs_close(hsq); 
}

int main(int argc, char **argv)
{
  lj_t *lj, *sq;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rc);
  sq = lj_open(N, D, rho, rc);

  lj_initsq(lj, ra, rb); /* we still need to pass ra & rb */
  lj->usesq = 0;
  lj_energy(lj);
  
  lj_initsq(sq, ra, rb);
  sq->esqinf = lj->esqinf = esqinf;
  
  domc(lj, sq);

  lj_close(lj);
  lj_close(sq);
  return 0;
}
