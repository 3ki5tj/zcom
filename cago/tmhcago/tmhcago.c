#define HAVE_REAL 1
/* typedef float real; */
typedef double real;

#include "time.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_CFG
#define ZCOM_AV
#define ZCOM_LOG
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
int nstcom = 10;

double tmh_ensexp = 0.0;
double tmh_emin, tmh_emax, tmh_de = 1.0;
int tmh_guesserange = 0; /* guess energy range */
double tmh_tps = 0.5; /* thermostat temperature */
double tmh_tp0 = 0.1, tmh_tp1 = 1.0, tmh_dtp = 0.0;
double tmh_erg0, tmh_erg1, tmh_derg = 1.0;
double tmh_ampmax = 1e-4, tmh_ampc = 2.0, tmh_lgvdt = 2e-3;

/* parameters for guess the erange */
int tmh_tctrun = -1, tmh_trep = 20000;

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
  CFGGETI(nstcom);

  CFGGETD(tmh_tp0);
  CFGGETD(tmh_tp1);
  CFGGETD(tmh_dtp);

  CFGGETD(tmh_de);
  CFGGETD(tmh_derg);

  CFGGETD(tmh_tps);

  /* if any of erg0, erg1, emin or emax is not given, guess energy range */
  CFGGETD(tmh_emin); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_emax); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_erg0); if (ret != 0) tmh_guesserange = 1;
  CFGGETD(tmh_erg1); if (ret != 0) tmh_guesserange = 1;
  CFGGETI(tmh_tctrun); CFGGETI(tmh_trep);

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

static int tmhrun(tmh_t *tmh, cago_t *go, double nsteps, double step0)
{
  int it = 0, stop = 0;
  double t, amp, ampadj, dhde;
  logfile_t *log = log_open("TRACE");

  ampadj = amp = tmh_ampmax;
  for (t = step0; t <= nsteps; t++) {
    dhde = tmh_getdhde(tmh, go->epot) * tmh_tps / tmh->tp;
    cago_vv(go, (real)dhde, mddt);
    if (it == 0) cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, (real)(tmh_tps), thermdt);
    go->rmsd = cago_rotfit(go, go->x, NULL);

    tmh_eadd(tmh, go->epot);
    tmh_dhdeupdate(tmh, go->epot, ampadj);

    if (++it % 10 == 0) {
      it = 0;
      tmh_tlgvmove(tmh, go->epot, tmh_lgvdt);
      /* update amplitude */
      if ((amp = tmh_ampc/t) > tmh_ampmax) amp = tmh_ampmax;
      ampadj = amp*tmh->tp0/tmh->tp;
    }
    if ((int) fmod(t, 10000) == 0) {
      printf("t %8g, T %6.3f, U %8.3f, rmsd %6.3f, dhde = %6.3f\n", t, tmh->tp, go->epot, go->rmsd, dhde);
    }

    if ((int)fmod(t, nsttrace) == 0) {
      log_printf(log, "%g %g %g %g %d %d %g\n",
          t, go->epot, tmh->tp, go->rmsd, tmh->iec, tmh->itp, dhde);
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
  log_close(log);
  return 0;
}

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(ARGOPT_LONGOPT); /* for -maxh */
  argopt_regarg(ao, NULL, &fncfg, "cfgfile");
  argopt_reghelp(ao, "-h");
  argopt_regversion(ao, "--version");
  argopt_regopt(ao, "-v", "%b", &verbose, "verbose");
  argopt_regopt(ao, "-maxh", "%lf", &maxtime, "max. simulation hours");
  argopt_parse(ao, argc, argv);
  if (argopt_set(ao, maxtime)) maxtime *= 3600.*.98;
  argopt_close(ao); 
}

int main(int argc, char **argv)
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
  printf("epotref %g\n", go->epotref);

  if (isctn) { /* load previous data */
    if (0 != cago_readpos(go, go->x, go->v, fnpos))
      return -1;
    if (0 != tmh_loaderange(fndhde, &tmh_tp0, &tmh_tp1, &dtp,
          &tmh_erg0, &tmh_erg1, &tmh_derg, &tmh_emin, &tmh_emax, &tmh_de,
          &tmh_ensexp, &tmh_dhdeorder))
      return -1;
  } else { /* fresh new run */
    if (tmh_guesserange) 
        printf("this feature is not supported\n"), exit(1);
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

