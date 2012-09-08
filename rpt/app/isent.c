/* computing a temperature profile for the 2D Ising model
 * using entropic sampling to cover the energy space */
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

double nsteps = 1e6; /* number of sweeps */
int nevery = 1;  /* compute temperatures every this number of steps */
char *fnprof = "profis.dat"; /* profile */
char *fnflow = "flowis.dat"; /* profile */

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation SWEEPS");
  argopt_add(ao, "-e", "%d", &nevery, "interval of computing temperatures");  
  argopt_add(ao, "-o", NULL, &fnprof, "name of the profile file (output)");
  argopt_add(ao, "-f", NULL, &fnflow, "name of the flow file (output)");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  nsteps *= N; /* convert to the number of steps */
  argopt_dump(ao);
  argopt_close(ao);
}

/* compute the reference beta */
static void getrefbet(double *refbet, const double *logdos)
{
  int i, len = L*L+1, il, ir;
  double *arr;
  const double ninf = -10000;
  
  /* smoothed DOS to remove ruggedness */
  xnew(arr, len);
  for (i = 0; i < len; i++) {
    if (logdos[i] > ninf) {
      arr[i] = logdos[i];
      continue;
    }
    for (il = i - 1; il >= 0; il--)
      if (logdos[il] > ninf)
        break;
    if (il < 0) { /* cannot smooth it */
      arr[i] = logdos[i];
      continue;
    }
    for (ir = i + 1; ir < len; ir++) 
      if (logdos[ir] > ninf)
        break;
    if (ir >= len) {
      arr[i] = logdos[i];
      continue;
    }
    /* interpolate */
    arr[i] = (logdos[il]*(ir-i) + logdos[ir]*(i-il))/(ir - il);
  }
  
  refbet[0] = (arr[1] - arr[0])/EDEL;
  refbet[len-1] = (arr[len-1] - arr[len-2])/EDEL;
  for (i = 1; i < len - 1; i++) {
    refbet[i] = (arr[i+1] - arr[i-1])/(2.*EDEL);
  }
  free(arr);
}

/* entropic sampling */
static int move(ising_t *is, int *de0)
{
  int id, h, acc = 0;
  int eo, en, ieo, ien;
  double dlogg;

  IS2_PICK(is, id, h); /* set id, and compute the # of opposite neighbors */
  eo = is->E; /* the current energy */
  en = is->E + h*2;  /* the new energy */
  *de0 = h * 2;
  if (h == 0) { /* a move that does not change energy */
    acc = 1;
  } else { /* judge from the DOS */
    ieo = (eo - EMIN)/EDEL;
    ien = (en - EMIN)/EDEL;     
    dlogg = is->logdos[ien] - is->logdos[ieo];
    /*printf("%d, %d, %g\n", eo, en, dlogg);*/
    if (dlogg <= 0 || rnd0() < exp(-dlogg))
      acc = 1;
  }
  if (acc) {
    IS2_FLIP(is, id, h);
    return 2*h;
  } else {
    return 0;
  }
}

/* compute the overall flow */
static int getflow(rpti_t **ta, double *logdos, const char *fn)
{
  double *f4, *f8, cnt, cntp, cntn, a, b;
  int i, j, ic = ta[0]->m/2;
  FILE *fp;
  
  xnew(f4, L*L+1);
  xnew(f8, L*L+1);
  for (i = 0; i < L*L; i++) {
    cnt = rpti_cnt(ta[i]);   if (cnt <= 0) cnt = 0;
    cntn = rpti_cnt(ta[i+1]); if (cntn <= 0) cntn = 0;
    a = ta[i]->h[ic+1]/cnt;
    if (a > 0) a = log(a); else a = -1000;
    a += logdos[i];
    b = ta[i+1]->h[ic-1]/cntn;
    if (b > 0) b = log(b); else b = -1000;
    b += logdos[i+1];
    f4[i] = lnadd(a, b);
    if (i > 0) {
      cntp = rpti_cnt(ta[i-1]); if (cntp <= 0) cntp = 0; 
      a = ta[i-1]->h[ic+2]/cntp;
      if (a > 0) a = log(a); else a = -1000;
      a += logdos[i-1];
      b = ta[i+1]->h[ic-2]/cntn;
      if (b > 0) b = log(b); else b = -1000;
      b += logdos[i+1];        
      f8[i] = lnadd(a, b);
    }
  }
  
  xfopen(fp, fn, "w", return -1);
  for (i = 1; i < L*L; i++) {
    fprintf(fp, "%d %g %g %g %g %g %g %g ", EMIN + i*EDEL,
      f4[i], f8[i], logdos[i-1], logdos[i], logdos[i+1],
      .5*(logdos[i]+logdos[i+1]), .5*(logdos[i-1]+logdos[i+1]) );
    cnt = rpti_cnt(ta[i]);
    if (cnt <= 0) cnt = 1;
    for (j = 0; j < ta[i]->m; j++)
      fprintf(fp, "%g ", ta[i]->h[j]/cnt);
    fprintf(fp, "\n");
  }
  fclose(fp);
  free(f4);
  free(f8);
  return 0;
}

static double *getserr(double *bet, double *ref, int n)
{
  int i;
  double *serr;
  
  xnew(serr, L*L+2);
  for (i = n/2; i <= n; i++)
    serr[i+1] = serr[i] + .5*EDEL*(bet[i] - ref[i] + bet[i+1] - ref[i+1]);  
  for (i = n/2; i > 0; i--)
    serr[i-1] = serr[i] - .5*EDEL*(bet[i] - ref[i] + bet[i-1] - ref[i-1]);
  return serr;
}

/* save the data */
static int savedata(int *hist, rpti_t **ta, double *refbet, const char *fn)
{
  FILE *fp;
  int i, j, n = L*L, verbose = 0;
  double bet0, bet1;
  double *beti, *bets0, *bets1, *beth0, *beth1, *serr, *serr0, *serr1, *serr2, *serr3;
  
  xfopen(fp, fn, "w", return -1);
  
  xnew(beti, L*L+1);
  xnew(bets0, L*L+1);
  xnew(bets1, L*L+1);
  xnew(beth0, L*L+1);
  xnew(beth1, L*L+1);  
  for (i = 0; i <= n; i++) {
    beti[i] = rpti_bet(ta[i]);
    bets0[i] = rpti_bets(ta[i], 0);
    bets1[i] = rpti_bets(ta[i], 1);
    beth0[i] = rpti_beth(ta[i], 0);
    beth1[i] = rpti_beth(ta[i], 1);    
  }
  serr  = getserr(beti,  refbet, n);
  serr0 = getserr(bets0, refbet, n);
  serr1 = getserr(bets1, refbet, n);
  serr2 = getserr(beth0, refbet, n);
  serr3 = getserr(beth1, refbet, n);
    
  for (i = 0; i <= n; i++) {
    fprintf(fp, "%d %d ", EMIN + i*EDEL, hist[i]);
    bet1 = rpti_bet1(ta[i], &bet0);
    fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f "
      "%.6f %.6f %.6f %.6f %.6f ", 
      bet0, bet1, beti[i], bets0[i], bets1[i], beth0[i], beth1[i], refbet[i],
      serr[i], serr0[i], serr1[i], serr2[i], serr3[i]);
    for (j = 0; j < ta[i]->m; j++)
      fprintf(fp, "%d ", ta[i]->h[j]);
    fprintf(fp, "\n");
    
    if (verbose >= 2) {  printf("E %d, bet %g, %g, %g;\n", 
    EMIN+i*EDEL, bet0, bet1, beti[i]); getchar(); }
  }
  free(beti);
  free(serr);
  fclose(fp);
  return 0;
}

/* run Monte Carlo simulation */
static int mc(ising_t *is, double nst)
{
  double st;
  int i, ie, de, de0;
  int *hist;
  rpti_t **ta;
  double *refbet;
  
  xnew(hist, L*L+1);
  xnew(ta, L*L+1);
  for (i = 0; i <= L*L; i++) ta[i] = rpti_open(-8, 8, 4, 0);
  
  ie = (is->E - EMIN)/EDEL;
  for (st = 0; st <= nst; st++) {
    hist[ie]++;
    
    de = move(is, &de0);
    
    /* no need to add the entry if de == 0 
       for the integral solution, but it is necessary
       for the formula 2de/de2 */
    rpti_add(ta[ie], de0);
    
    if (de != 0) ie += de/4;
  }
  //printf("simlation finished; press enter.\n"); getchar();
  
  xnew(refbet, L*L+1);
  getrefbet(refbet, is->logdos);
  savedata(hist, ta, refbet, fnprof);
  getflow(ta, is->logdos, fnflow);
  
  /* free data */
  for (i = 0; i <= L*L; i++) rpti_close(ta[i]);
  free(ta);
  free(hist);
  free(refbet);
  return 0;
}

int main(int argc, char **argv)
{
  ising_t *is;

  doargs(argc, argv);
  die_if ((is = is2_open(L)) == NULL, "cannot open is2, L %d\n", L);

  /* read in the density of states */
  die_if (is2_loadlogdos(is, NULL) != 0,
      "no dos for the %dx%d model\n", L, L);

  mc(is, nsteps);
      
  is2_close(is);
  return 0;
}

