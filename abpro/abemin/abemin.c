/* check energy and force */
#define HAVE_REAL 1
//typedef float real;
typedef double real;

#define ZCOM_PICK
#define ZCOM_ABPRO
#include "zcom.h"

const char *fn = NULL;
int itmax = 1000000;
double tol = 1e-10;
int overwrite = 0;
int verbose = 0;

static int getinfo(const char *fn, int *d, int *model, int *seqid)
{
  FILE *fp;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  if (3 != fscanf(fp, "# %d %d %d", d, model, seqid)) {
    fprintf(stderr, "cannot scan basic info from %s\n", fn);
    fclose(fp);
    return -1;
  }
  fclose(fp);
  return 0;
}

/* print help and die */
static void help(const char *prog)
{
  printf("%s [OPTIONS] file\n", prog);
  printf("check energy of a position file\n\n");
  printf("OPTIONS:\n");
  printf(" -n: followed by max. # of iterations, def %d\n", itmax);
  printf(" -t: followed by tolerance, def %g\n", tol);
  printf(" -w: overwrite the file, if significantly minimized\n");
  exit(1);
}

/* handle arguments */
static int doargs(int argc, const char **argv)
{
  int i, ch;
  const char *val;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fn = argv[i];
      continue;
    }
    ch = argv[i][1];
    if (strchr("nt", ch) != NULL) {
      if (i == argc - 1) help(argv[0]);
      val = argv[++i];
    }

    if (ch == 'n') {
      itmax = atoi(val);
    } else if (ch == 't') {
      tol = atof(val);
    }
    if (ch == 'w') overwrite = 1;
    else if (ch == 'v') verbose = 1;
  }
  if (fn == NULL) help(argv[0]);
  return 0;
}
  
int main(int argc, const char **argv)
{
  abpro_t *ab;
  int seqid, d, model, flags;
  real Em, E;

  doargs(argc, argv);
  if (getinfo(fn, &d, &model, &seqid) != 0) 
    return -1;
  ab = ab_open(seqid, d, model, 0.);
  if (verbose) fprintf(stderr, "load %s\n", fn);
  if (ab_readpos(ab, ab->x, NULL, fn) != 0) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  E = ab_energy(ab, ab->x, 0);
  flags = AB_LMREGISTER | (verbose ? AB_VERBOSE : 0);
  Em = ab_localmin(ab, ab->x, itmax, tol, flags);
  fprintf(stderr, "E: %.8f -> %.8f\n", E, Em);
  if (fabs(E - Em) > 1e-6 && overwrite) {
    printf("update energy...\n");
    ab_writepos(ab, ab->lmx, NULL, fn);
  }
  ab_close(ab);
  return 0;
}

