/* check energy and force */
#define HAVE_REAL 1
typedef double real;

#define ZCOM_PICK
#define ZCOM_ABPRO
#define ZCOM_ARGOPT
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
  argopt_regopt(ao, "-t", "%f", &tol, "tolerance of energy minimizer");
  argopt_regopt(ao, "-N", "%d", &sh_itmax, "max. # of iterations of SHAKE");
  argopt_regopt(ao, "-T", "%f", &sh_tol, "SHAKE tolerance");
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
    ab_initconstr(ab, localevel);
    ab_updconstr(ab);
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

