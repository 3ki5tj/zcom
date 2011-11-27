#define HAVE_REAL 1
/* typedef float real; */
typedef double real;

#include "time.h"
#define ZCOM_PICK
#define ZCOM_CFG
#define ZCOM_ARGOPT
#define ZCOM_LOG
#define ZCOM_ABPRO
#define ZCOM_AV
#define ZCOM_TMH
#include "zcom.h"

const char *fntp = "tmhab.t", *fndhde = "tmhab.e", *fnehis = "tmhab.ehis";
const char *fnpos = "ab.pos";
const char *fncfg = "tmhab.cfg";
const char *fnlog = "tmhab.tr";
double nsteps = 1000000*1000;
int seqid = 10, d = 3, model = 2;
int tmh_dhdeorder = 1;
double tmh_dhdemin = 0.1, tmh_dhdemax = 10.0;
real mddt = 5e-3f;
real thermdt = 5e-3f;
int usebrownian = 0;
real brdt = 5e-4;
int lgclevel = 0; /* local geometry optimization: aba */
unsigned milcshake = AB_MILCSHAKE;

double tmh_ensexp = 0.0;
double tmh_emin, tmh_emax, tmh_de = 1.0;
int tmh_guesserange = 0; /* guess energy range */
double tmh_tps = 0.3; /* thermostat temperature */
double tmh_tp0 = 0.1, tmh_tp1 = 1.0, tmh_dtp = 0.0;
double tmh_erg0, tmh_erg1, tmh_derg = 1.0;
double tmh_elimit = 1e9;
double tmh_springk = 0.0;
double tmh_ampmax = 1e-4, tmh_ampc = 2.0, tmh_ampmin = 0.0;
double tmh_entampmax = 1e-2, tmh_entampc = 100.0, tmh_entampmin = 0.0;
int tmh_updampc = 0;
double tmh_lgvdt = 2e-3;

/* parameters for guess the erange */
int tmh_annealcnt = 1;
int tmh_teql = 400000, tmh_tctrun = 400000, tmh_trep = 20000;
double tmh_erg0margin = 0., tmh_erg1margin = 0.1;
int tmh_entropic = 0; /* entropic sampling */

double nsttrace = 1e4;
double nstsave = 1e6;
int nstmv = 10;

double maxtime = 1e9;
time_t time0, time1;

int tmh_srand = 0; /* seed for initial random configuration */

int verbose = 0;

#define CFGGETI(x) { cfg_add(cfg, #x, "%d",  &x, #x); }
#define CFGGETD(x) { cfg_add(cfg, #x, "%lf", &x, #x); }
#define CFGGETR(x) { cfg_add(cfg, #x, "%r",  &x, #x); }

/* load setting from configuration */
static int loadcfg(const char *fn)
{
  cfg_t *cfg;

  die_if ((cfg = cfgopen(fn)) == NULL, "cannot open %s\n", fn);

  cfg_add(cfg, "seqid", "!%d", &seqid, "sequence id");
  cfg_add(cfg, "d", "!%d", &d, "dimension");
  cfg_add(cfg, "model", "!%d", &model, "model (1 or 2)");

  CFGGETD(nsteps);
  CFGGETR(mddt);
  CFGGETR(thermdt);
  CFGGETI(usebrownian);
  CFGGETR(brdt);
  CFGGETI(lgclevel);

  CFGGETD(tmh_tp0);
  CFGGETD(tmh_tp1);
  CFGGETD(tmh_dtp);

  CFGGETD(tmh_de);
  CFGGETD(tmh_derg);

  CFGGETD(tmh_tps);
  
  /* if any of erg0, erg1, emin or emax is not given, guess energy range */
  CFGGETD(tmh_emin);
  CFGGETD(tmh_emax);
  CFGGETD(tmh_erg0);
  CFGGETD(tmh_erg1);
  CFGGETI(tmh_annealcnt);
  CFGGETD(tmh_erg0margin);
  CFGGETD(tmh_erg1margin);
  CFGGETI(tmh_teql); CFGGETI(tmh_tctrun); CFGGETI(tmh_trep);

  CFGGETD(tmh_elimit);
  CFGGETD(tmh_springk);

  CFGGETI(tmh_dhdeorder);
  CFGGETD(tmh_dhdemin);
  CFGGETD(tmh_dhdemax);
  CFGGETD(tmh_ensexp);
  CFGGETD(tmh_ampmax);
  CFGGETD(tmh_ampmin);
  CFGGETD(tmh_ampc);
  CFGGETI(tmh_updampc);
  CFGGETD(tmh_lgvdt);
  
  CFGGETI(tmh_entropic);
  CFGGETD(tmh_entampmax);
  CFGGETD(tmh_entampmin);
  CFGGETD(tmh_entampc);

  CFGGETI(tmh_srand);

  CFGGETD(nsttrace);
  CFGGETD(nstsave);
  cfg_add(cfg, "nstmv", "%d", &nstmv, "number of steps per tempering");

  die_if (cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE) != 0, "failed match %d\n", fn);
  if (!(cfg_set(cfg, tmh_emin) && cfg_set(cfg, tmh_emax) 
        && cfg_set(cfg, tmh_erg0) && cfg_set(cfg, tmh_erg1)))
    tmh_guesserange = 1;
  if (tmh_srand) srand((unsigned) tmh_srand);
  cfg_close(cfg);
  return 0;
}

/* run constant temperature run */
static double ctrun(abpro_t *ab, double tp, double *edev,
    int teql, int tmax, int trep)
{
  int t;
  av_t av[1];
  double eav, dev;

  ab->t = 0.; /* reset time */
  av_clear(av);
  for (t = 1; t <= teql+tmax; t++) {
    if (usebrownian) {
      ab_brownian(ab, (real)tp, 1.f, (real)brdt, AB_SOFTFORCE|milcshake);
    } else {
      ab_vv(ab, (real)( tmh_tps/tp ), (real)mddt, AB_SOFTFORCE|milcshake);
      if (t % 10 == 0) ab_rmcom(ab, ab->x, ab->v);
      ab_vrescale(ab, (real) tmh_tps, (real)thermdt);
    }
    if (t % 1000 == 0 && ab->lgcon && ab->lgact < ab->lgcnt)
      ab_updconstr(ab, 0);
    if (t > teql) av_add(av, ab->epot);
    if (trep > 0 && t % trep == 0) {
      fprintf(stderr, "t = %g, tp = %g, epot %g, ekin %g, lgc %d/%d\n", 
          ab->t, tp, ab->epot, ab->ekin, ab->lgact, ab->lgcnt);
    }
  }
  eav = av_getave(av);
  dev = av_getdev(av);
  printf("tp %g, eav %g, dev %g\n", tp, eav, dev);
  if (edev != NULL) *edev = dev;
  return eav;
}

/* tmh run */
static int tmhrun(tmh_t *tmh, abpro_t *ab, double nsteps, double step0)
{
  int it = 0, stop = 0;
  double t, dhde;
  logfile_t *log = log_open(fnlog);

  for (t = step0; t <= nsteps; t++) {
    dhde = tmh_getdhde(tmh, ab->epot) * tmh_tps / tmh->tp;
    if (usebrownian == 2) {
      ab_brownian(ab, (real) (tmh_tps * dhde), 1, (real) brdt, AB_SOFTFORCE|milcshake);
    } else if (usebrownian == 1) {
      ab_brownian(ab, (real) tmh_tps, (real) dhde, (real) brdt, AB_SOFTFORCE|milcshake);
    } else {
      ab_vv(ab, (real) dhde, (real) mddt, AB_SOFTFORCE|milcshake);
      if (it == 0) ab_rmcom(ab, ab->x, ab->v);
      ab_vrescale(ab, (real) tmh_tps, (real) thermdt);
    }
    
    if (ab->epot < ab->emin + 0.05 || (tmh->iec < 3 && rnd0() < 1e-4)) {
      double em = ab->emin;
      int lgcon = ab->lgcon;
      ab->lgcon = 0; /* turn off LGC during energy minimization */
      if (ab_localmin(ab, ab->x, 0, 0., 0, 0., AB_LMREGISTER|AB_LMWRITE) < em)
        printf("emin = %10.6f from %10.6f t %g tp %g.%30s\n", ab->emin, ab->epot, t, tmh->tp, "");
      ab->lgcon = lgcon;
    }
    if (++it % nstmv == 0) {
      it = 0;
      if (ab->lgcon && ab->lgact < ab->lgcnt)
        ab_updconstr(ab, 0);
      if (tmh_entropic) {
        tmh_ezmoves(tmh, ab->epot, 1.0);
      } else { /* tweak updating factor by tps/tp */
        tmh_ezmove(tmh, ab->epot, tmh_tps/tmh->tp, tmh_lgvdt);
      }
    }
  
    if ((int) fmod(t, nsttrace) == 0) {
      log_printf(log, "%g %g %g %d %d %g\n",
          t, ab->epot, tmh->tp, tmh->iec, tmh->itp, dhde);
      if (verbose) fprintf(stderr, "t = %g epot = %g, tp = %g, iec %d, itp %d, dhde %g;%20s\r",
          t, ab->epot, tmh->tp, tmh->iec, tmh->itp, dhde, "");
      time(&time1);
      if (difftime(time1, time0) > maxtime) {
        printf("\n\nsimulation exceeds %g hours\n\n", maxtime/3600.);
        stop = 1;
      }
    }
    if ((int) fmod(t, nstsave) == 0 || stop) {
      mtsave(NULL);
      tmh_save(tmh, fntp, fnehis, fndhde, tmh->wl->lnf, t);
      ab_writepos(ab, ab->x, ab->v, fnpos);
    }
    if (stop) break;
  }
  log_close(log);
  return 0;
}

/* guess energy range from constant temperature simulation */
static void guess_erange(abpro_t *ab, double tp0, double tp1, 
    double *erg0, double *erg1, double *emin, double *emax,
    int annealcnt, int teql, int trun, int trep)
{
  double x, edv0, edv1;
  int i;

  *erg1 = ctrun(ab, tp1, &edv1, teql, trun, trep);
  for (i = 1; i <= annealcnt; i++) {
    x = 1.*i/(annealcnt + 1);
    ctrun(ab, x*tp0 + (1. - x)*tp1, NULL, teql, trun, trep);
  }
  *erg0 = ctrun(ab, tp0, &edv0, teql, trun, trep);
  ab_localmin(ab, ab->x, 0, 0., 0, 0., AB_LMREGISTER|AB_LMWRITE);

  x = *erg1 - *erg0;
  *emin = *erg0 - x;
  *emax = *erg1 + 2*x;
  *erg0 += x*tmh_erg0margin;
  *erg1 -= x*tmh_erg1margin;
  printf("Guessing energy range...\n%g: %g%+g; %g: %g%+g E (%g, %g), emin = %g\n", 
      tp0, *erg0, edv0, tp1, *erg1, edv1, *emin, *emax, ab->emin); 
}

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  double maxh = 1e9;
  argopt_t *ao = argopt_open(ARGOPT_LONGOPT); /* for -h */
  argopt_regopt(ao, "-maxh", "%lf", &maxh, "maximal hours");
  argopt_regopt(ao, "-cfg", NULL, &fncfg, "configuration file");
  argopt_regopt(ao, "-g", NULL, &fnlog, "log file");
  argopt_reghelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  maxtime = maxh * 3600 * 0.95;
  argopt_close(ao);
}

int main(int argc, char **argv)
{
  abpro_t *ab;
  tmh_t *tmh;
  double step0 = 0, amp, dtp;
  int isctn; 

  doargs(argc, argv);
  time(&time0);
  die_if (loadcfg(fncfg) != 0, "cannot load %s\n", fncfg);
  isctn = fexists(fndhde); /* continue if fndhde */
  ab = ab_open(seqid, d, model, 0.1);
  if (lgclevel) {
    ab_initconstr(ab, lgclevel);
    milcshake = 0;
  }
  if (isctn) { /* load previous data */
    if (0 != ab_readpos(ab, ab->x, ab->v, fnpos))
      return -1;
    if (0 != tmh_loaderange(fndhde, &tmh_tp0, &tmh_tp1, &dtp, 
          &tmh_erg0, &tmh_erg1, &tmh_derg, &tmh_emin, &tmh_emax, &tmh_de,
          &tmh_ensexp, &tmh_dhdeorder))
      return -1;
  } else { /* fresh new run */
    if (tmh_guesserange)
      guess_erange(ab, tmh_tp0, tmh_tp1, &tmh_erg0, &tmh_erg1, &tmh_emin, &tmh_emax,
          tmh_annealcnt, tmh_teql, tmh_tctrun, tmh_trep);
    /* equilibrate to thermostat temperature */
    ctrun(ab, tmh_tps, NULL, tmh_teql, tmh_tctrun, tmh_trep);
    ab->epot = ab_energy(ab, ab->x, 0);
    printf("finish equilibration, epot %g, ekin %g\n", ab->epot, ab->ekin);
  }

  tmh = tmh_open(tmh_tp0, tmh_tp1, tmh_dtp, tmh_erg0, tmh_erg1, tmh_derg, 
      tmh_emin, tmh_emax, tmh_de, tmh_ensexp, tmh_dhdeorder);
  tmh->elimit = tmh_elimit;
  tmh->springk = tmh_springk;
  if (tmh_entropic) {
    tmh_initwlcvg(tmh, tmh_entampc, tmh_entampmax, sqrt(0.1), 0.0, 
      tmh_entampmin, tmh_updampc ? WLCVG_UPDLNFC : 0);
  } else {
    tmh_initwlcvg(tmh, tmh_ampc, tmh_ampmax, sqrt(0.1), 0.0, 
      tmh_ampmin, tmh_updampc ? WLCVG_UPDLNFC : 0);
  }
  printf("erange (%g, %g), active (%g, %g)\n", tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);
  tmh->dhdemin = tmh_dhdemin;
  tmh->dhdemax = tmh_dhdemax;

  tmh_settp(tmh, tmh_tps); /* tp = tps */
  die_if (isctn && tmh_load(tmh, fnehis, fndhde, &amp, &step0) != 0,
    "cannot load tmh, %s, %s\n", fndhde, fnehis);

  tmhrun(tmh, ab, nsteps, step0);

  tmh_close(tmh);
  ab_close(ab);
  return 0;
}

