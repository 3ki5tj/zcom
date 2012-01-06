/* merge data from three data sets */
#define ZCOM_PICK
#define ZCOM_DISTR
#define ZCOM_CFG
#include "zcom.h"

#define S 3

const char *fncfg = "ene.cfg";

const char *fnds  = "data/dsmerge.dat";
const char *fndsarr[S] = {"datal/ds_c.dat", "data/ds_c.dat", "datah/ds_c.dat"};

const char *fndsb  = "data/dsbmerge.dat";
const char *fndsbarr[S] = {"datal/dsb_c.dat", "data/dsb_c.dat", "datah/dsb_c.dat"};

const char *fnfe = "fe.dat"; /* free energy */

double dbet[S] = {1/.8 - 1, 0, 1/1.2 - 1};
double wt[S] = {0.01, 1.0, 0.5};
/* we assign the Tl with a smaller weight as it has a larger correlation time */

double wtb[S] = {1.0, 1.0, 1.0};

double emin = -1400, emax = -1200, edel = 0.1;

int iitype = 0; /* 0: Adib-Jarsynski; 1: modulated */
int halfwin = 50; /* half window size */
int halfwinb = 50; /* half window size for db */
int mfhalfwin = 0; /* half window size for mean force */

/* adaptive window parameters */
double gam = 1.0;
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
  cfg_add(cfg, "halfwinb", "%d", &halfwinb, "half of the number of bins in each side of the window, "
      "for the fractional identity; 0: guess, -1: adaptive");
  cfg_add(cfg, "iitype", "%d", &iitype, "integral identity type: 0: Adib-Jarzynski, 1: modulated");
  cfg_add(cfg, "mfhalfwin", "%d", &mfhalfwin, "half of the number of bins in the window, "
      "for mean force; 0: single bin, < 0: plain average");
  cfg_add(cfg, "gamma", "%lf", &gam, "half window amplification factor");  
  cfg_add(cfg, "mlimit", "%d", &mlimit, "maximal # of bins in adaptive window");  
  cfg_add(cfg, "sampmin", "%lf", &sampmin, "minimal number of samples to estimate sig(mf)");
  cfg_match(cfg, 0);
  cfg_close(cfg);
}

/* for canonical ensemble only
 * merging and compute mean force */
INLINE void mfmhc(distr_t *d, distr_t **darr, double *dbeta, double *wtk, int narr)
{
  int i, m = 0, k, n = d->n;
  double s, sf, sf2, sdv;

  for (i = 0; i < n; i++) {
    distrsum_t *ds;
    s = sf = sf2 = sdv = 0;
    for (k = 0; k < narr; k++) {
      ds = darr[k]->arr + i;
      s += wtk[k] * ds->s;
      sf += wtk[k] * (ds->sf + ds->s * dbeta[k]);
      sf2 += wtk[k] * (ds->sf2 + 2 * ds->sf * dbeta[k] + dbeta[k] * dbeta[k] * ds->s);
      sdv += wtk[k] * (ds->sdv);
    }
    ds = d->arr + i;
    ds->s = s;
    ds->sf = sf;
    ds->sf2 = sf2;
    ds->sdv = sdv;
    /* extend to both sides, stop as soon as we have one data point */
    for (m = 1; m < n && s <= 0.; m++) {
      if (i - m >= 0) { /* extend to the left */
        for (k = 0; k < narr; k++) {
          ds = darr[k]->arr + i - m;
          s += wtk[k] * ds->s;
          sf += wtk[k] * ds->sf + ds->s * dbeta[k];
          sdv += wtk[k] * ds->sdv;
        }
      }
      if (i + m < n) { /* extend to the right */
        for (k = 0; k < narr; k++) {
          ds = darr[k]->arr + i + m;
          s += wtk[k] * ds->s;
          sf += wtk[k] * ds->sf + ds->s * dbeta[k];
          sdv += wtk[k] * ds->sdv;
        }
      }
    }
    d->mf[i] = (s > 0.) ? (sf / s) : 0.;
    d->mdv[i] = (s > 0.) ? (sdv / s) : 0.;
  }
}

/* estimate lnz from d->rho */
INLINE void estlnz(distr_t *d, double *dbeta, double *dlnz, int narr)
{
  int i, k, n = d->n;
  double x, dx = d->dx;

  for (k = 0; k < narr; k++) {
    dlnz[k] = -10000.0;
    for (i = 0; i <= n; i++) {
      x = d->xmin + i * dx;
      if (d->rho[i] <= 0.) continue; /* no contribution */
      dlnz[k] = lnadd(dlnz[k], -dbeta[k] * x + log(d->rho[i] * dx));
    }
  }
}

INLINE void iimhc(distr_t *d, distr_t **darr, double *dbeta, double *wt,
    double *dlnz, int narr)
{
  int i, j, k, jl, jr, n = d->n, iter;
  double x, den, denk, *totk, tot, tot1, dx = d->dx;

  /* construct ln(rho) from the mean force */
  distr_iimf(d);

  xnew(totk, narr);

  for (d->hsum[0] = 0., i = 0; i < n; i++) {
    for (tot = 0, k = 0; k < narr; k++) {
      tot += x = wt[k] * darr[k]->arr[i].s;
      totk[k] += x;
    }
    d->hsum[i + 1] = d->hsum[i] + tot;
  }
  tot = d->hsum[n]; /* total number of visits */
  
  for (iter = 0; iter < 50; iter++) {
    estlnz(d, dbeta, dlnz, narr);
    tot1 = 0.;
    /* estimate using integral identity */
    for (i = 0; i <= n; i++) {
      /* compute # of visits in a window around i */
      jl = d->jl[i];
      jr = d->jr[i];
      d->his[i] = d->hsum[jr] - d->hsum[jl];

      for (den = 0, k = 0; k < narr; k++) {
        /* integrate over exp(lnrho[j] - lnrho[i])
         * Note: we can loop from j0 to j1, with little difference */
        for (denk = 0., j = jl; j <= jr; j++) {
          x = exp(d->lnrho[j] - d->lnrho[i] - dbeta[k] * (j - i) * dx);
          /* if (d->norm) x *= d->norm[j]; not RDF, no need */
          denk += (j == jl || j == jr) ? x * .5 : x;
        }
        denk *= totk[k] * dx * exp(-dbeta[k] * (d->xmin + i * dx) - dlnz[k]);
        den += denk;
      }
      d->rho[i] = (den > 0.) ? d->his[i] / den : 0.;
      d->his[i] /= tot * dx * (jr - jl);
      tot1 += d->rho[i] * dx;
    }
   
    for (i = 0; i <= n; i++) d->rho[i] /= tot1;
  }
  free(totk);
  printf("tot1 %g; ", tot1);
  for (k = 0; k < S; k++) printf("dlnz[%d] %g, ", k, dlnz[k]);
  printf("\n");
}

/* reweight histogram, mean force to make it look like they are from a single temperature simulation */
INLINE void decorate(distr_t *d, distr_t **darr, double *dbeta, double *wtk,
    double *dlnz, int narr)
{
  int i, m = 0, k, n = d->n;
  double s, sf, sf2, sdv, wt;

  for (i = 0; i < n; i++) {
    distrsum_t *ds;
    s = sf = sf2 = sdv = 0;
    for (k = 0; k < narr; k++) {
      ds = darr[k]->arr + i;
      s += wtk[k] * ds->s;
      sf += wtk[k] * (ds->sf + ds->s * dbeta[k]);
      sf2 += wtk[k] * (ds->sf2 + 2 * ds->sf * dbeta[k] + dbeta[k] * dbeta[k] * ds->s);
      sdv += wtk[k] * ds->sdv;
    }
    for (wt = 0, k = 0; k < narr; k++) {
      wt += wtk[k] * exp(-dbeta[k]*(d->xmin + (i + .5) * d->dx) - dlnz[k]);
    }
    ds = d->arr + i;
    ds->s = s / wt;
    ds->sf = sf / wt;
    ds->sf2 = sf2 / wt;
    ds->sdv = sdv / wt;
  }
}


/* merge distribution in *darr into one */
static void merge(distr_t *d, distr_t **darr, double *dbeta, double *wt,
    double *dlnz, int narr, int m)
{
  int i, s, n = d->n;

  mfmhc(d, darr, dbeta, wt, narr); /* computing the mean force */
  if (m < 0) m = 0;
  distr_winfixed(d, gam, m);
  iimhc(d, darr, dbeta, wt, dlnz, narr);
  decorate(d, darr, dbeta, wt, dlnz, narr);
}

/* free energy */
INLINE void fene(distr_t *d, distr_t *db,
    double tp0, double tp1, double tp, const char *fn)
{
  double bet, bet0 = 1.0/tp1, bet1 = 1./tp0, dbet = 0.01;
  int n = (int)((bet1 - bet0)/dbet + .5), i;
  double *barr, *dlnz, *dlnzb;
  FILE *fp;
  
  xnew(barr, n + 1);
  xnew(dlnz, n + 1);
  xnew(dlnzb, n + 1);
  for (i = 0; i < n; i++)
    barr[i] = bet0 + dbet * i - 1.0/tp;
  estlnz(d, barr, dlnz, n);
  estlnz(db, barr, dlnzb, n);
  xfopen(fp, fn, "w", exit(1));
  for (i = 0; i < n; i++) {
    bet = barr[i] + 1.0/tp;
    fprintf(fp, "%g %g %g\n", 1.0/bet, -dlnzb[i]/bet, -dlnz[i]/bet);
  }
  fclose(fp);
  free(barr);
  free(dlnz);
  free(dlnzb);
}

int main(void)
{
  int s;
  distr_t *darr[S], *d, *dbarr[S], *db;
  double dlnz[S];

  loadcfg(fncfg);
  d  = distr_open(emin, emax, edel);
  for (s = 0; s < S; s++) {
    darr[s] = distr_open(emin, emax, edel);
    die_if(0 != distr_load(darr[s], fndsarr[s]), 
        "failed to load data from %s\n", fndsarr[s]);
  }
  
  merge(d, darr, dbet, wtb, dlnz, S, halfwin);

  distr_save(d, fnds);

  /* repeat the same thing on dsb */
  db = distr_open(emin, emax, edel);
  for (s = 0; s < S; s++) {
    dbarr[s] = distr_open(emin, emax, edel);
    die_if(0 != distr_load(dbarr[s], fndsbarr[s]), 
        "failed to load data from %s\n", fndsbarr[s]);
  }
  
  merge(db, dbarr, dbet, wtb, dlnz, S, halfwinb);

  distr_save(db, fndsb);
  fene(d, db, 0.7, 1.5, 1.0, fnfe);

  /* close everything */
  distr_close(d);
  for (s = 0; s < S; s++)
    distr_close(darr[s]);
  
  distr_close(db);
  for (s = 0; s < S; s++)
    distr_close(dbarr[s]);

  return 0;
}

