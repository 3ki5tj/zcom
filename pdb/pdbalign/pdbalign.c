/* align backbone structures */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ROTFIT
#define ZCOM_PDB
#include "zcom.h"

const char *prog = "pdbalign";
const char *fn1 = NULL;
const char *fn2 = NULL;
const char *fnout = "fit.pdb";
enum {ALN_CA = 0, ALN_BB = 1, ALN_LAST};
int aligntype = 1;
int verbose = 1;


/* print help and die */
static void help(void)
{
  printf("%s [OPTIONS] ref.pdb x.pdb\n", prog);
  printf("Rotate & translate x.pdb to a best fitting postion to ref.pdb\n\n");
  printf("OPTIONS:\n");
  printf(" -o: followed by the output file, default \"%s\"\n", fnout);
  printf(" -a: followed by alignment type 0: C-alpha, 1: backbone, default %d\n", aligntype);
  printf(" -v: followed by the verbose level, default %d\n", verbose);
  printf(" -h: print this message\n\n");
  exit(1);
}

/* handle arguments */
static int doargs(int argc, const char **argv)
{
  int i, j, ch, acnt = 0;
  const char *val;

  prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      if (acnt == 0) {
        fn1 = argv[i];
      } else if (acnt == 1) {
        fn2 = argv[i];
      } else help();
      acnt++;
      continue;
    }
    ch = argv[i][1];
    if (strchr("aov", ch)) { /* argument options */
      val = argv[i] + 2;
      if (*val == '\0') {
        if (i == argc - 1) {
          fprintf(stderr, "need arg. after %s\n", argv[i]);
          help();
        }
        val = argv[++i];
      }
      if (ch == 'o') {
        fnout = val;
      } else if (ch == 'a') {
        aligntype = atoi(val);
        if (aligntype < 0 || aligntype >= ALN_LAST) help();
      } else if (ch == 'v') {
        verbose = atoi(val);
      } else {
        fprintf(stderr, "program error: -%c is not handled\n", ch);
      }
    } else { /* simple options */
      for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
        if (ch == 'h') help();
      }
    }
  }
  if (fn1 == NULL || fn2 == NULL) help();
  printf("fitting %s against %s, type: %d\n", fn2, fn1, aligntype);
  return 0;
}

static rv3_t *load(const char *fn, int *pn, pdbmodel_t **pmdl)
{
  pdbmodel_t *m;
  pdbaac_t *c;
  int i, n;
  rv3_t *x;

  die_if ((m = pdbm_read(fn, verbose)) == NULL,
      "cannot read pdb %s\n", fn);
  die_if ((c = pdbaac_parse(m, verbose)) == NULL,
      "failed to parse the pdb %s\n", fn);
  *pn = n = c->nres;
  if (aligntype == ALN_CA) {
    xnew(x, n);
    for (i = 0; i < n; i++)
     rv3_copy(x[i],   c->res[i].xca); 
  } else if (aligntype == ALN_BB) {
    xnew(x, 3*n);
    for (i = 0; i < n; i++) {
     rv3_copy(x[3*i],   c->res[i].xn); 
     rv3_copy(x[3*i+1], c->res[i].xca); 
     rv3_copy(x[3*i+2], c->res[i].xc); 
    }
  } else fatal("unknown alignment type %d\n", aligntype);
  *pmdl = m;
  pdbaac_free(c);
  return x;
}

static void pdbm_rottrans(pdbmodel_t *m, real rot[3][3], real trans[3])
{
  real x1[3];
  int i;

  for (i = 0; i < m->n; i++) {
    rv3_rot(x1, rot, m->at[i].x);
    rv3_inc(x1, trans);
    rv3_copy(m->at[i].x, x1);
  }
}

/* assign weights (mass) for each atom */
static real *mkweights(int alntype, int nr)
{
  int i, n;
  real *w = NULL;
  
  if (alntype == ALN_CA) return NULL;
  
  if (alntype == ALN_BB) {
    n = nr*3;
    xnew(w, n);
    for (i = 0; i < nr; i++) {
      w[3*i] = 14; /* N */
      w[3*i+1] = 12; /* CA */
      w[3*i+2] = 12; /* C */
    }
  }
  return w;
}

int main(int argc, const char **argv)
{
  pdbmodel_t *m1, *m2;
  real *w;
  rv3_t *x1, *x2;
  real rot[3][3], trans[3];
  int nr, n1, n2, n;

  doargs(argc, argv);

  x1 = load(fn1, &n1, &m1);
  x2 = load(fn2, &n2, &m2);
  if (n1 != n2) {
    fprintf(stderr, "Warning: %d != %d\n", n1, n2);
  }
  nr = (n1 < n2) ? n1 : n2;
  if (aligntype == ALN_CA) {
    n = nr;
  } else if (aligntype == ALN_BB) {
    n = nr*3;
  } 
  w = mkweights(aligntype, nr);
  rotfit3(x2, NULL, x1, w, n, rot, trans);
  if (verbose) {
    mat3_print(rot, "Rotation", "%8.3f", 1);
    rv3_print(trans, "Translation", "%8.3f", 1);
  }
  pdbm_rottrans(m2, rot, trans);
  pdbm_write(m2, fnout);

  pdbm_free(m1);
  pdbm_free(m2);
  if (w) free(w);
  free(x1);  free(x2);
  return 0;
}

