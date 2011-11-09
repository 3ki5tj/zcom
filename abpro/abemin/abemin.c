/* check energy and force */
#define HAVE_REAL 1
typedef double real;

#define ZCOM_PICK
#define ZCOM_ABPRO
#define ZCOM_ARGOPT
#define ZCOM_AV
#include "zcom.h"

const char *fnpos = NULL;
int itmax = 1000000;
double tol = 1e-10;
int overwrite = 0;
int verbose = 0;
int milcshake = 0;
int sh_itmax = 100000;
double sh_tol = 1e-10;
int localgeom = 0;
int localevel = 1;

/* get dimension and model from file */
static int getinfo(const char *fn, int *d, int *model, int *seqid)
{
  FILE *fp;

  xfopen(fp, fn, "r", return -1);
  if (3 != fscanf(fp, "# %d %d %d", d, model, seqid)) {
    fprintf(stderr, "cannot scan basic info from %s\n", fn);
    fclose(fp);
    return -1;
  }
  fclose(fp);
  return 0;
}

static void doargs(int argc, char **argv)
{
  int Overwrite = 0;
  argopt_t *ao = argopt_open(0);
  ao->desc = "Energy minimize a structure";
  argopt_regarg(ao, NULL, &fnpos, "input.pos");
  argopt_regopt(ao, "-n", "%d", &itmax, "max. # of iterations");
  argopt_regopt(ao, "-t", "%lf", &tol, "tolerance of energy minimizer");
  argopt_regopt(ao, "-N", "%d", &sh_itmax, "max. # of iterations of SHAKE");
  argopt_regopt(ao, "-T", "%lf", &sh_tol, "SHAKE tolerance");
  argopt_regopt(ao, "-w", "%b", &overwrite, "overwrite the file, if significantly minimized");
  argopt_regopt(ao, "-W", "%b", &Overwrite, "always overwrite the file");
  argopt_regopt(ao, "-m", "%b", &milcshake, "use MILCSHAKE");
  argopt_regopt(ao, "-q", "%d", &localevel, "local geometry level");
  argopt_regopt(ao, "-l", "%b", &localgeom, "display local geometry info.");
  argopt_regopt(ao, "-v", "%b", &verbose, "be verbose");
  argopt_reghelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  if (Overwrite) overwrite = 2;
  argopt_close(ao);
}

/* print contact */
static void ab_printcontact(abpro_t *ab)
{
  int i, di, n = ab->n, d = ab->d, cnt = 0;
  real dr, dr0 = (real) pow(2.0, 1.0/6);
  av_t av[1];

  for (di = 2; di <= 3; di++) { /* (i, i + di) */
    av_clear(av);
    if (di == 2) printf("ABA:\n"); else printf("ABBA\n");
    for (i = 0; i < n - di; i++) {
      if (ab->type[i] != 0 || ab->type[i + di] != 0) continue;
      if (d == 3) dr = rv3_dist(ab->x + i*3, ab->x + (i + di)*3);
      else dr = rv2_dist(ab->x + i*2, ab->x + (i + di)*2);
      av_add(av, dr);
      printf("%3d: %3d - %3d: %9.6f - %9.6f = %9.6f\n", ++cnt, i, i + di, dr, dr0, dr - dr0);
    }
    printf("average distance %9.6f +/- %9.6f\n", 
        av_getave(av), av_getdev(av));
  }
}

int main(int argc, char **argv)
{
  abpro_t *ab;
  int seqid, d, model, flags;
  real Em, E;

  doargs(argc, argv);
  if (getinfo(fnpos, &d, &model, &seqid) != 0) 
    return -1;
  ab = ab_open(seqid, d, model, 0.);
  if (verbose) fprintf(stderr, "load %s\n", fnpos);
  if (ab_readpos(ab, ab->x, NULL, fnpos) != 0) {
    fprintf(stderr, "cannot open %s\n", fnpos);
    return -1;
  }
  memcpy(ab->x1, ab->x, ab->n*d*sizeof(real));
  if (milcshake) 
    ab_milcshake(ab, ab->x1, ab->x, NULL, 0., sh_itmax, sh_tol, 0);
  else
    ab_shake(ab, ab->x1, ab->x, NULL, 0., sh_itmax, sh_tol, 0);

  E = ab_energy(ab, ab->x, 0);
  flags = AB_LMREGISTER | (verbose ? AB_VERBOSE : 0) | (milcshake ? AB_MILCSHAKE : 0);
  Em = ab_localmin(ab, ab->x, itmax, tol, sh_itmax, sh_tol, flags);
  fprintf(stderr, "E: %.8f -> %.8f\n", E, Em);
  if ((fabs(E - Em) > 1e-6 && overwrite) || overwrite > 1) {
    printf("update configuration...\n");
    ab_writepos(ab, ab->lmx, NULL, fnpos);
  }
  if (localgeom) {
    ab_printcontact(ab);
    ab_initconstr(ab, localevel); /* initialize constraints */
    ab_updconstr(ab); /* add constraints */
    printf("established %d/%d constraints, model %d, %dD\n", ab->lgact, ab->lgcnt, ab->model, ab->d);
    flags = verbose ? AB_VERBOSE : 0;
    E = ab_localmin(ab, ab->x, 1000, tol, 100, sh_tol, flags);
    ab->lgcon = 0; /* turn off */
    flags |= milcshake ? AB_MILCSHAKE : 0;
    Em = ab_localmin(ab, ab->x, itmax, tol, sh_itmax, sh_tol, flags);
    printf("local geometry: E %.8f (with constraints) -> %.8f\n", E, Em);
  }
  ab_close(ab);
  return 0;
}

