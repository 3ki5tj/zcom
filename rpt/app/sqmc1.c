/* random perturbation temperature profile with Monte Carlo simulation
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
int iumin = -1000;
int iumax = 0;
int iucnt; 
int nequil = 1000;
int nsteps = 1000000;

int dumax = 1000; /* histogram dimension */
char *fnbetp = "betp.dat"; /* energy-increment distribution */

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
  argopt_add(ao, "-m", "%r", &amp,    "amplitude of a MC move"); 
  argopt_add(ao, "-[", "%d", &iumin,  "minimal potential energy");
  argopt_add(ao, "-]", "%d", &iumax,  "maximal potential energy");
  argopt_add(ao, "-U", "%d", &dumax,  "maximal energy change for histogram");
  argopt_add(ao, "-d", NULL, &fnbetp, "output file for the beta profile");
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
static void domc(lj_t *lj, rpti_t **rpt)
{
  int t, acc = 0;
  int iu, idu;
  double u, du;
  static av_t avU[1], avu[1];

  for (t = 0; t < nequil; t++) /* warm up */
    lj_metro3d(lj, amp, 1.0f/tp);

  for (t = 1; t <= nsteps; t++) { /* real simulation */
    iu = lj->iepot;
    acc += lj_metrosq3d(lj, amp, 1.0f/tp);
    idu = lj->iepot - iu;
    rpti_add(rpt[iu - iumin], idu);
    av_add(avu, (real) idu);
    av_add(avU, lj->epot);
    //if (acc) {   printf("%d, %d, %g", idu, lj->iepot, lj->epot); getchar(); }
  }
  u = av_getave(avU)/lj->n;
  du = av_getdev(avu);

  printf("epot %g, acc %g, dudev %g\n",
     u, 1.*acc/nsteps, du);
}

static int savedata(rpti_t **rpt, const char *fn)
{
  double *his, *bp, *bp0, *bp1, *bpi;
  int i, i0 = -1, i1 = -1, ie;
  FILE *fp;
  
  xnew(his, iucnt + 1);
  xnew(bp0, iucnt + 1);
  xnew(bp1, iucnt + 1);
  xnew(bpi, iucnt + 1);
  xnew(bp,  iucnt + 1);  
  xfopen(fp, fn, "w", return -1);
  for (i = 0; i < iucnt; i++) {
    his[i] = rpti_cnt(rpt[i]);
    if (his[i] > 1e-14) { /* skip empty ones */
      if (i0 < 0) i0 = i;
      if (i > i1) i1 = i;
      bp1[i] = rpti_bet1(rpt[i], &bp0[i]);
      bpi[i] = rpti_bet(rpt[i]);
      if (fabs(bp0[i]) > .99*RPT_INF) bp0[i] = 0.f;
      if (fabs(bp1[i]) > .99*RPT_INF) bp1[i] = 0.f;
      if (fabs(bpi[i]) > .99*RPT_INF) bpi[i] = 0.f;
    }
  }

  /* differentiate the histogram */
  for (i = i0; i < i1; i++) {
    if (his[i] <= 0. || his[i+1] <= 0.) {
      bp[i] = 0.;
    } else {
      bp[i] = log(his[i+1]/his[i]);
    }
  }

  /* output */
  for (i = i0; i <= i1; i++) {
    ie = iumin + i;
    fprintf(fp, "%4d %.6f %.6f %.6f %.6f %.1f %.6f\n",
      ie, bp0[i], bp1[i], bpi[i], bp[i], ie + .5, his[i]);
  }
  fclose(fp);
  free(his);
  free(bp0);
  free(bp1);
  free(bpi);
  free(bp);

  return 0;
}

int main(int argc, char **argv)
{
  lj_t *lj;
  rpti_t **rpt;
  int i;

  doargs(argc, argv);
  lj = lj_open(N, D, rho, rb);
  lj_initsq(lj, ra, rb);
  xnew(rpt, iumax - iumin + 1);
  iucnt = iumax - iumin + 1;
  for (i = 0; i < iucnt; i++)
    rpt[i] = rpti_open(-dumax, dumax, 1, RPTI_HASINF);
   
  domc(lj, rpt);
  
  savedata(rpt, fnbetp);
  //rpti_wdist(rpt, fnehis);
  for (i = 0; i < iucnt; i++)
    rpti_close(rpt[i]);
  free(rpt);
  lj_close(lj);
  return 0;
}
