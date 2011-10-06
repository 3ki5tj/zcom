#define ZCOM_PICK
#define ZCOM_SS
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

static void help(const char *prog)
{
  printf("%s [OPTIONS] [proj]\n", prog);
  printf("histogram reweighting for tmh\n\n");
  printf("OPTIONS:\n");
  printf(" -p: followed by project name, def = %s\n", proj);
  printf(" -e: followed by dhde file\n");
  printf(" -h: followed by histogram file\n");
  printf(" -z: followed by lnz (partition function) file\n");
  printf(" -g: followed by lng (density of states) file\n");
  printf(" -n: followed by maximal number of iterations, def = %d\n", itmax);
  printf(" -t: followed by the tolerance of lnz convergence, def = %g\n", tol);
  printf(" -h: print this message\n");  
  exit(1);
}

static int doargs(int argc, char **argv)
{
  int i, j, ch;
  const char *val;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      proj = argv[i];
      continue;
    }
    ch = argv[i][1];
    if (strchr("pehzgnt", ch) != NULL) { /* argument options */
      val = argv[i] + 2;
      if (*val == '\0') {
        if (i == argc - 1) {
          fprintf(stderr, "need arg. after %s\n", argv[i]);
          help(argv[0]);
        }
        val = argv[++i];
      }
      if (ch == 'p') {
        proj = val;
      } else if (ch == 'e') {
        fndhde = ssdup(val);
      } else if (ch == 'h') {
        fnehis = ssdup(val);
      } else if (ch == 'z') {
        fnlnz = ssdup(val);
      } else if (ch == 'g') {
        fndos = ssdup(val);
      } else if (ch == 'n') {
        itmax = atoi(val);
      } else if (ch == 't') {
        tol = atof(val);
      } else {
        fprintf(stderr, "program error: -%c is not handled\n", ch);
      }
    } else { /* simple options */
      for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
        if (ch == 'v') verbose = 1;
        else if (ch == 'h') help(argv[0]);
      }
    }
  }
  if (fndhde == NULL) { fndhde = ssdup(proj); sscat(fndhde, ".e"); }
  if (fnehis == NULL) { fnehis = ssdup(proj); sscat(fnehis, ".ehis"); }
  if (fnlnz  == NULL) { fnlnz = ssdup(proj);  sscat(fnlnz,  ".lnz"); }
  if (fndos  == NULL) { fndos = ssdup(proj);  sscat(fndos,  ".lng"); }
  return 0;
}

int main(int argc, char **argv)
{
  tmh_t *tmh;
  double tp0, tp1, dtp, emin, emax, de, erg0, erg1, derg, ensexp = 2.0;
  int dhdeorder = 1;
  double amp, t0;

  doargs(argc, argv);

  if (0 != tmh_loaderange(fndhde, &tp0, &tp1, &dtp, &erg0, &erg1, &derg,
        &emin, &emax, &de, &ensexp, &dhdeorder))
    return -1;

  tmh = tmh_open(tp0, tp1, dtp, erg0, erg1, derg, emin, emax, de,
      ensexp, dhdeorder);

  if (tmh_load(tmh, fnehis, fndhde, &amp, &t0) != 0) {
    fprintf(stderr, "cannot load tmh\n");
    return -1;
  }
  tmh_calcdos(tmh, itmax, tol, fndos, fnlnz);
  tmh_close(tmh);

  return 0;
}
