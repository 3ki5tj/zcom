#define HAVE_REAL 1
/* typedef float real; */
typedef double real;

#include "time.h"
#define ZCOM_PICK
#define ZCOM_CFG
#define ZCOM_AV
#define ZCOM_TRACE
#define ZCOM_CAGO
#define ZCOM_TMH
#include "zcom.h"

/* ca-go parameters */
real cago_kb = 200.f; /* bond */
real cago_ka = 40.f;  /* angle */
real cago_kd1 = 1.f;  /* dihedral phi*1 */
real cago_kd3 = .5f;  /* dihedral phi*3 */
real cago_nbe = 1.f;
real cago_nbc = 4.f;
real cago_rc = 6.f; /* cutoff distance */

const char *fntp = "tmhcago.t", *fndhde = "tmhcago.e", *fnehis = "tmhcago.ehis";
char *fnpdb = NULL; /* initial structure */
const char *fnpos = "cago.pos"; /* position file */
const char *fncfg = "tmhcago.cfg";

double nsteps = 1000000*100;
int tmh_dhdeorder = 1;
double tmh_dhdemin = 0.1, tmh_dhdemax = 10.0;
real mddt = 2e-3f;
real thermdt = 2e-3f;

double tmh_ensexp = 0.0;
double tmh_emin, tmh_emax, tmh_de = 1.0;
int tmh_guesserange = 0; /* guess energy range */
double tmh_tps = 0.5; /* thermostat temperature */
double tmh_rmsd0 = 2.0, tmh_rmsd1 = 20.0, tmh_drmsd = 0.5;
double tmh_tp0 = 0.1, tmh_tp1 = 1.0, tmh_dtp = 0.0;
double tmh_erg0, tmh_erg1, tmh_derg = 1.0;
double tmh_ampmax = 1e-4, tmh_ampc = 2.0, tmh_lgvdt = 2e-3;

/* parameters for guess the erange */
int tmh_annealcnt = 1;
int tmh_teql = 400000, tmh_tctrun = 400000, tmh_trep = 20000;
double tmh_erg0margin = 0., tmh_erg1margin = 0.1;

double nsttrace = 1e4;
double nstsave = 1e6;

double maxtime = 1e9;
time_t time0, time1;

int tmh_srand = 0; /* seed for initial random configuration */

const char *prog = "tmhgo";
int verbose = 0;

#define CFGGETI(x) { ret = cfgget(cfg, &x, #x, "%d");  printf("%s: %d\n", #x, x); }
#define CFGGETD(x) { ret = cfgget(cfg, &x, #x, "%lf"); printf("%s: %g\n", #x, x); }
#define CFGGETR(x) { ret = cfgget(cfg, &x, #x, rfmt);  printf("%s: %g\n", #x, x); }

/* load setting from configuration */
static int loadcfg(const char *fn)
{
  cfgdata_t *cfg;
  int ret;
  const char *rfmt; /* format string for real */

  die_if ((cfg = cfgopen(fn)) == NULL, "cannot open %s\n", fn);

  rfmt = ((sizeof(real) == sizeof(double)) ? "%lf" : "%f");

  ret = cfgget(cfg, &fnpdb, "pdb", "%s");
  printf("pdb file [%s]\n", fnpdb);

  /* ca-go parameters */
  CFGGETR(cago_kb);
  CFGGETR(cago_ka);
  CFGGETR(cago_kd1);
  CFGGETR(cago_kd3);
  CFGGETR(cago_nbe);
  CFGGETR(cago_nbc);
  CFGGETR(cago_rc);

  CFGGETD(nsteps);

  CFGGETR(mddt);
  CFGGETR(thermdt);

  CFGGETD(tmh_rmsd0);
  CFGGETD(tmh_rmsd1);
  CFGGETD(tmh_drmsd);
  //CFGGETD(tmh_tp0);
  //CFGGETD(tmh_tp1);
  //CFGGETD(tmh_dtp);

  CFGGETD(tmh_de);
  CFGGETD(tmh_derg);

  CFGGETD(tmh_tps);

  /* if any of erg0, erg1, emin or emax is not given, guess energy range */
  CFGGETD(tmh_emin); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_emax); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_erg0); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_erg1); if (ret != 0) tmh_guesserange = 1;
  CFGGETI(tmh_annealcnt);
  CFGGETD(tmh_erg0margin);
  CFGGETD(tmh_erg1margin);
  CFGGETI(tmh_teql); CFGGETI(tmh_tctrun); CFGGETI(tmh_trep);

  CFGGETI(tmh_dhdeorder);
  CFGGETD(tmh_dhdemin);
  CFGGETD(tmh_dhdemax);
  CFGGETD(tmh_ensexp);
  CFGGETD(tmh_ampmax);
  CFGGETD(tmh_ampc);
  CFGGETD(tmh_lgvdt);

  CFGGETI(tmh_srand);
  if (tmh_srand) srand((unsigned) tmh_srand);

  CFGGETD(nsttrace);
  CFGGETD(nstsave);

  cfgclose(cfg);
  return 0;
}

/* run constant temperature run */
static double ctrun(cago_t *go, double tp, av_t *vrmsd, av_t *vepot,
    int teql, int tmax, int trep)
{
  int t;

  av_clear(vepot);
  av_clear(vrmsd);
  go->t = 0.; /* reset time */
  for (t = 1; t <= teql+tmax; t++) {
    cago_vv(go, 1.f, mddt);
    if (t % 10 == 0) cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, (real) tp, thermdt);
    go->rmsd = cago_rotfit(go, go->x, NULL);
    if (t > teql) {
      av_add(vepot, go->epot);
      av_add(vrmsd, go->rmsd);
    }
    if (trep > 0 && t % trep == 0) {
      fprintf(stderr, "t = %g, tp = %g, epot %g, ekin %g\n", go->t, tp, go->epot, go->ekin);
    }
  }
  printf("tp %g, epot %g(%g), rmsd %g(%g)\n", tp, 
      av_getave(vepot), av_getdev(vepot), av_getave(vrmsd), av_getdev(vrmsd));
  return av_getave(vepot);
}

static int tmhrun(tmh_t *tmh, cago_t *go, double nsteps, double step0)
{
  int it = 0, stop = 0;
  double t, amp, dhde;

  amp = tmh_ampmax;
  for (t = step0; t <= nsteps; t++) {
    dhde = tmh_getdhde(tmh, tmh->ec, tmh->iec)*tmh_tps/tmh->tp;
    cago_vv(go, (real)dhde, mddt);
    if (it == 0) cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, (real)(tmh_tps), thermdt);

    tmh_eadd(tmh, go->epot);
    tmh_dhdeupdate(tmh, go->epot, amp);

    if (++it % 10 == 0) {
      it = 0;
      tmh_tlgvmove(tmh, go->epot, tmh_lgvdt);
      /* update amplitude */
      if ((amp = tmh_ampc/t) > tmh_ampmax) amp = tmh_ampmax;
    }

    if ((int)fmod(t, nsttrace) == 0) {
      wtrace_buf("%g %g %g %d %d %g\n",
          t, go->epot, tmh->tp, tmh->iec, tmh->itp, dhde);
      if (verbose) fprintf(stderr, "t = %g epot = %g, tp = %g, iec %d, itp %d, dhde %g;%20s\r",
          t, go->epot, tmh->tp, tmh->iec, tmh->itp, dhde, "");
      time(&time1);
      if (difftime(time1, time0) > maxtime) {
        printf("\n\nsimulation exceeds %g hours\n\n", maxtime/3600.);
        stop = 1;
      }
    }
    if ((int)fmod(t, nstsave) == 0 || stop) {
      mtsave(NULL);
      tmh_save(tmh, fntp, fnehis, fndhde, amp, t);
      cago_writepos(go, go->x, go->v, fnpos);
    }
    if (stop) break;
  }
  /* finish */
  wtrace_buf(NULL, NULL);
  return 0;
}

/* guess energy range from constant temperature simulation */
static void guess_erange(cago_t *go,
    double rmsd0, double rmsd1,
    double *tp0, double *tp1,
    double *erg0, double *erg1, double *emin, double *emax,
    int annealcnt, int teql, int trun, int trep)
{
  double x, edv0, edv1;
  int i;
  av_t vep[1], vrmsd[1], vtp[1];

  cago_cvgmdrun(go, rmsd0, 10, 0.01, sqrt(0.1), 0.1, 
      vtp, vep, 1.0, 0.01, 100.0, teql, trep);
  exit(1);
/*
  *erg1 = ctrun(go, tp0, rmsd0, &edv1, teql, trun, trep);
  for (i = 1; i <= annealcnt; i++) {
    x = 1.*i/(annealcnt + 1);
    ctrun(go, x*tp1 + (1. - x)*tp0, NULL, teql, trun, trep);
  }
  *erg0 = ctrun(go, tp0, &edv0, teql, trun, trep);

  x = *erg1 - *erg0;
  *emin = *erg0 - x;
  *emax = *erg1 + 2*x;
  *erg0 += x*tmh_erg0margin;
  *erg1 -= x*tmh_erg1margin;
  printf("Guessing energy range...\n%g: %g%+g; %g: %g%+g E (%g, %g)\n",
      tp0, *erg0, edv0, tp1, *erg1, edv1, *emin, *emax);
*/
}

/* print usage and die */
static void help(void)
{
  printf("%s [OPTIONS] [your.cfg]\n", prog);
  printf("OPTIONS:\n"
      "  -maxh: followed by maximal hour\n");
  exit(1);
}

/* handle input arguments */
static int doargs(int argc, const char **argv)
{
  int i;

  prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fncfg = argv[i];
      continue;
    }
    if (strcmp(argv[i], "-maxh") == 0) {
      if (i == argc - 1) help();
      i++;
      maxtime = atof(argv[i])*3600.*0.98;
    } else if (strcmp(argv[i], "-v") == 0) {
      verbose = strlen(argv[i] + 1); /* number of v's, e.g., -vvv --> verbose = 3 */
    } else if (strcmp(argv[i], "-h") == 0) {
      help();
    } else {
      printf("unknown option %s\n", argv[i]);
      help();
    }
  }
  return 0;
}

int main(int argc, const char **argv)
{
  cago_t *go;
  tmh_t *tmh;
  double step0 = 0, amp, dtp;
  int isctn;

  doargs(argc, argv);
  time(&time0);

  die_if (loadcfg(fncfg) != 0, "cannot load %s\n", fncfg);

  /* assume a continuation run if fndhde exists */
  isctn = fexists(fndhde);

  if ((go = cago_open(fnpdb, cago_kb, cago_ka, cago_kd1, cago_kd3,
          cago_nbe, cago_nbc, cago_rc)) == NULL) {
    fprintf(stderr, "error initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, thermdt); /* also init a MD system */

  if (isctn) { /* load previous data */
    if (0 != cago_readpos(go, go->x, go->v, fnpos))
      return -1;
    if (0 != tmh_loaderange(fndhde, &tmh_tp0, &tmh_tp1, &dtp,
          &tmh_erg0, &tmh_erg1, &tmh_derg, &tmh_emin, &tmh_emax, &tmh_de,
          &tmh_ensexp, &tmh_dhdeorder))
      return -1;
  } else { /* fresh new run */
    if (tmh_guesserange) {
      guess_erange(go,
          tmh_rmsd0, tmh_rmsd1,
          &tmh_tp0, &tmh_tp1,
          &tmh_erg0, &tmh_erg1,
          &tmh_emin, &tmh_emax,
          tmh_annealcnt, tmh_teql, tmh_tctrun, tmh_trep);
    }
    go->epot = cago_force(go, go->f, go->x);
    printf("finish equilibration, epot %g, ekin %g\n", go->epot, go->ekin);
  }

  tmh = tmh_open(tmh_tp0, tmh_tp1, tmh_dtp, tmh_erg0, tmh_erg1, tmh_derg,
      tmh_emin, tmh_emax, tmh_de, tmh_ensexp, tmh_dhdeorder);
  printf("erange (%g, %g), active (%g, %g)\n",
      tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);
  tmh->dhdemin = tmh_dhdemin;
  tmh->dhdemax = tmh_dhdemax;

  if (isctn && tmh_load(tmh, fnehis, fndhde, &amp, &step0) != 0) {
    fprintf(stderr, "cannot load tmh\n");
    return -1;
  } else {
    tmh_settp(tmh, tmh_tps);
  }

  tmhrun(tmh, go, nsteps, step0);

  tmh_close(tmh);
  cago_close(go);
  return 0;
}

