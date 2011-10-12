#define ZCOM_PICK
#define ZCOM_SS
#define ZCOM_ARGOPT
#define ZCOM_TMH
#include "zcom.h"

const char *proj = "tmhpt";
char *fndhde = NULL;
char *fnehis = NULL;
char *fndos = NULL;
char *fnlnz = NULL;
int verbose = 0;
int itmax = 10000;
double tol = 1e-6;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "Histogram reweighting for tmh";
  argopt_regarg(ao, NULL, &proj, "project");
  argopt_regopt(ao, "-p", NULL, &proj, "project");
  argopt_regopt(ao, "-e", NULL, &fndhde, "dhde file");
  argopt_regopt(ao, "-h", NULL, &fnehis, "histogram file");
  argopt_regopt(ao, "-z", NULL, &fnlnz, "lnz (partition function) file");
  argopt_regopt(ao, "-g", NULL, &fndos, "lng (density of states) file");
  argopt_regopt(ao, "-n", "%d", &itmax, "maximal number of iterations");
  argopt_regopt(ao, "-t", "%lf", &tol, "tolerance of lnz convergence");
  argopt_reghelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
  if (fndhde == NULL) { fndhde = ssdup(proj); sscat(fndhde, ".e"); }
  if (fnehis == NULL) { fnehis = ssdup(proj); sscat(fnehis, ".ehis"); }
  if (fnlnz  == NULL) { fnlnz = ssdup(proj);  sscat(fnlnz,  ".lnz"); }
  if (fndos  == NULL) { fndos = ssdup(proj);  sscat(fndos,  ".lng"); }
}

int main(int argc, char **argv)
{
  tmh_t *tmh;
  double tp0, tp1, dtp, emin, emax, de, erg0, erg1, derg, ensexp = 2.0;
  int dhdeorder = 1;
  double amp, t0;

  doargs(argc, argv);

  if (0 != tmh_loaderange(fndhde, &tp0, &tp1, &dtp, &erg0, &erg1, &derg,
        &emin, &emax, &de, &ensexp, &dhdeorder)) {
    ssdelall();
    return -1;
  }

  tmh = tmh_open(tp0, tp1, dtp, erg0, erg1, derg, emin, emax, de,
      ensexp, dhdeorder);

  if (tmh_load(tmh, fnehis, fndhde, &amp, &t0) != 0) {
    fprintf(stderr, "cannot load tmh\n");
    tmh_close(tmh);
    ssdelall();
    return -1;
  }
  tmh_calcdos(tmh, itmax, tol, fndos, fnlnz);
  tmh_close(tmh);
  ssdelall(); /* free all strings */
  return 0;
}
