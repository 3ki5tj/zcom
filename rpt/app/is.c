/* computing a temperature profile for the 2D Ising model
 * constant temperature simulation */
#define  IS2_LB  5 /* L == 2^LB */

/* local parameters for the 2D Ising model */
#define L (1 << IS2_LB)
#define N (L*L)
#define EMIN (-2*N)
#define EMAX (2*N)
#define EDEL 4

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ISING2
#define ZCOM_RPT
#define ZCOM_ARGOPT
#include "zcom.h"

double tp = 1.0f;
int nequil = 10000;
double nsteps = 1e5*N;
int nevery = 1;  /* compute temperatures every this number of steps */
char *fnehis = "ehismc.dat"; /* energy-increment distribution */
char *fnehisd = NULL; /* adjusted energy-increment distribution */

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-T", "%lf", &tp,    "temperature");
  argopt_add(ao, "-0", "%d", &nequil, "number of equilibration");
  argopt_add(ao, "-1", "%d", &nsteps, "number of simulation sweeps");
  argopt_add(ao, "-e", "%d", &nevery, "interval of computing temperatures");  
  argopt_add(ao, "-o", NULL, &fnehis, "output file for the energy-increment histogram");
  argopt_add(ao, "-O", NULL, &fnehisd,"output file for the adjusted energy-increment histogram");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

/* Metropolis */
static int move(ising_t *is, int *de0)
{
  int id, h, acc = 0;

  IS2_PICK(is, id, h); /* set id, and compute the # of opposite neighbors */
  *de0 = h * 2;
  if (h <= 0 || mtrand() < is->uproba[h])
  //if (metroacc1(2*h, 1.0f/tp))
    acc = 1;
 
  if (acc) {
    IS2_FLIP(is, id, h);
    return 2*h;
  } else {
    return 0;
  }
}

/* run Monte Carlo simulation */
static int mc(ising_t *is, double bet, rpti_t *rpt, rpti_t *rptd)
{
  int i, de, de0;
  double U, du, bp0, bp1, bpi, bps0, bps1, bph0, bph1, bpd0, bpd1, bpdi;
  static av_t avU[1], avu[1];  
  
  IS2_SETPROBA(is, bet);
  for (i = 0; i <= nequil; i++)
    move(is, &de0);
    
  for (i = 1; i <= nsteps; i++) {
    de = move(is, &de0);
    
    rpti_add(rpt, de0);
    rpti_add(rptd, de);
    av_add(avu, de);
    av_add(avU, is->E);
  }
  U = av_getave(avU);
  du = av_getdev(avu);
  bp1 = rpti_bet1(rpt, &bp0);
  bpi = rpti_bet(rpt);
  bps0 = rpti_bets(rpt, 0);
  bps1 = rpti_bets(rpt, 1);
  bph0 = rpti_beth(rpt, 0);
  bph1 = rpti_beth(rpt, 1);
  bpd1 = rpti_bet1(rptd, &bpd0);
  bpdi = rpti_bet(rptd);
  printf("epot %g, du %g, bp0 %.6f, bp1 %.6f, bpi %.6f, bps0 %.6f, bps1 %.6f, "
    "bph0 %.6f, bph1 %.6f, bpd0 %.6f, bpd1 %.6f, bpdi %.6f\n",
     U, du, bp0, bp1, bpi, bps0, bps1, 
     bph0, bph1, bpd0, bpd1, bpdi);
  return 0;
}

int main(int argc, char **argv)
{
  ising_t *is;
  rpti_t *rpt, *rptd;

  doargs(argc, argv);
  if ((is = is2_open(L)) == NULL) {
    fprintf(stderr, "cannot init is2\n");
    return -1;
  }
  rpt = rpti_open(-8, 8, 4, 0);
  rptd = rpti_open(-8, 8, 4, 0);

  mc(is, 1.0f/tp, rpt, rptd);

  rpti_wdist(rpt, fnehis);
  if (fnehisd) rpti_wdist(rptd, fnehisd);
  rpti_close(rpt);
  rpti_close(rptd);      
  is2_close(is);
  return 0;
}
