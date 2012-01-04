/* load two databases, 1: sample, 2: reference,
 * apply different windows to the sample
 * compute error w.r.t. the reference */
#define ZCOM_PICK
#define ZCOM_DISTR
#define ZCOM_CFG
#include "zcom.h"

const char *fncfg = "ene.cfg";
const char *fnds = "data/dsmerge.dat"; /* reference */
//const char *fnds = "data/ds_c_ii.dat"; /* reference */
const char *fndsb = "data/dsb_c_ii.dat"; /* sample */

double emin = -1380, emax = -1240, edel = 0.1;

int iitype = 0; /* 0: Adib-Jarsynski; 1: modulated */
int halfwin = 50; /* half window size */
int mfhalfwin = 0; /* half window size for mean force */

/* adaptive window parameters */
int mlimit = -1; /* maximal # of bins for a window in adaptive window */
double sampmin = 400; /* minimal number of samples to estimate sig(mf) */

/* load parameters from .cfg file */
static void loadcfg(const char *fncfg)
{
  cfg_t *cfg = cfg_open(fncfg);
  if (cfg == NULL) return;
  cfg_add(cfg, "emin", "%lf", &emin, "minimal energy");
  cfg_add(cfg, "emax", "%lf", &emax, "maximal energy");
  cfg_add(cfg, "edel", "%lf", &edel, "energy interval");
  cfg_add(cfg, "halfwin", "%d", &halfwin, "half of the number of bins in each side of the window, "
      "for the fractional identity; 0: guess, -1: adaptive");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, "
      "for mean force; 0: single bin, < 0: plain average");
  cfg_add(cfg, "mlimit", "%d", &mlimit, "maximal # of bins in adaptive window");  
  cfg_add(cfg, "sampmin", "%lf", &sampmin, "minimal number of samples to estimate sig(mf)");
  cfg_match(cfg, 0);
  cfg_close(cfg);
}

/* compute error */
static void calcerr(const distr_t *db, const distr_t *d,
    double *hre, double *hrp, double *hrs, double *rre, double *rrp, double *rrs)
{
  int i, n = d->n;
  double *his0, *hisref, smpsiz;
  
  xnew(hisref, n + 1);
  xnew(his0, n + 1)
  for (smpsiz = 0, i = 0; i < n; i++) {
    hisref[i] = d->arr[i].s;
    his0[i] = db->arr[i].s;
    smpsiz += db->arr[i].s;
  }

  /* error of the histogram */
  *hre = ksdif1d(his0, d->rho, n, DSDIF_HIS1, smpsiz, hrp, NULL); 
  *hrs = entdif1d(his0, d->rho, n, DSDIF_HIS1, NULL); 
  /* error of the distribution */
  *rre = ksdif1d(db->rho, d->rho, n, 0, smpsiz, rrp, NULL);
  *rrs = entdif1d(d->rho, db->rho, n, DSDIF_HIS1, NULL); 
  free(his0);
  free(hisref);
}

/* scan window size and report error */
static void winscan(distr_t *db, distr_t *d)
{
  int m;
  double hre, hrp, hrs, rre, rrp, rrs, are, arp, ars;

  /* adaptive window */
  /* distr_iiez(db, 1, -1, 0, mlimit, sampmin); */
/*
  m = 0;
  distr_iiez(db, 1, m, 0, mlimit, sampmin);
  calcerr(db, d, &hre, &hrp, &hrs, &rre, &rrp, &rrs);
  printf("%d %g %g %g %g %g %g\n", m,
      hre, hrp, hrs, rre, rrp, rrs);
*/

  for (m = 5; m <= 250; m += 5) {
    distr_iiez(db, 0, m, 0, mlimit, sampmin);
    calcerr(db, d, &hre, &hrp, &hrs, &are, &arp, &ars);
    distr_iiez(db, 1, m, 0, mlimit, sampmin);
    calcerr(db, d, &hre, &hrp, &hrs, &rre, &rrp, &rrs);
    printf("%d %g %g %g %g %g %g %g %g %g\n", m,
        hre, hrp, hrs, rre, rrp, rrs, are, arp, ars);
  }
}

int main(void)
{
  distr_t *d, *db;

  loadcfg(fncfg);
  d = distr_open(emin, emax, edel);
  db = distr_open(emin, emax, edel);
  die_if(0 != distr_load(d, fnds), "failed to load data from %s\n", fnds);
  die_if(0 != distr_load(db, fndsb), "failed to load data from %s\n", fndsb);
  
  /* construct a best guess */
  if (strstr(fnds, "merge") == 0) /* only for single distribution */
    distr_iiez(d, 1, halfwin, mfhalfwin, mlimit, sampmin);
  winscan(db, d);

  distr_close(d);
  distr_close(db);
  return 0;
}

