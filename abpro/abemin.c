/* check energy and force */
#define HAVE_REAL 1
//typedef float real;
typedef double real;
#include "abpro.c"

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

int main(int argc, const char **argv)
{
  abpro_t *ab;
  int seqid, d, model;
  const char *fn;
  real Em, E;

  if (argc < 2) {
    fprintf(stderr, "%s file\n", argv[0]);
    return -1;
  }
  fn = argv[1];
  if (getinfo(fn, &d, &model, &seqid) != 0) 
    return -1;
  ab = ab_open(seqid, d, model, 0.);
  fprintf(stderr, "load %s\n", fn);
  if (ab_readpos(ab, ab->x, NULL, fn) != 0) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  E = ab_energy(ab, ab->x, 0);
  Em = ab_localmin(ab, ab->x, 1000000, 1e-12);
  fprintf(stderr, "E: %.8f -> %.8f\n", E, Em);
  if (fabs(E - Em) > 1e-6) {
    printf("update energy...\n");
    ab_writepos(ab, ab->lmx, NULL, fn);
  }
  ab_close(ab);
  return 0;
}

