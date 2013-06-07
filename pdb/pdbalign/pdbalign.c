/* align backbone structures */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ROTFIT
#define ZCOM_PDB
#define ZCOM_ARGOPT
#include "zcom.h"

const char *fn1 = NULL;
const char *fn2 = NULL;
const char *fnout = "fit.pdb";
enum {ALN_CA = 0, ALN_BB = 1, ALN_MC = 2, ALN_HEAVY = 8, ALN_ALL = 9, ALN_LAST = 10};
int aligntype = 1;
int verbose = 1;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "Align structures in file1 and file2 by rotating and translating that in file2";
  argopt_regarg(ao, "!", &fn1, "file1");
  argopt_regarg(ao, "!", &fn2, "file2");
  argopt_regopt(ao, "-o", NULL, &fnout, "fit2");
  argopt_regopt(ao, "-a", "%d", &aligntype, "type 0: C-alpha, 1: backbone, 2: mainchain, 8: heavy, 9: all");
  argopt_regopt(ao, "-v", "%d", &verbose, "verbose level");
  argopt_reghelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  if (aligntype < 0 || aligntype >= ALN_LAST) argopt_help(ao);
  argopt_close(ao);
}

/* only simple elements */
static real elem2mass(const char *elem)
{
  static struct { const char *name; double mass; } elemtab[] = {
    {"H", 1.00794}, {"D", 2.0141}, {"C", 12.0107}, {"N", 14.0067},
    {"O", 15.9994}, {"F", 18.9984}, {"Na", 22.98976928}, {"Mg", 24.3050f},
    {"Al", 26.9815386}, {"Si", 28.0855}, {"P", 30.973762}, {"S", 32.065},
    {"Cl", 35.453}, {"K", 39.0983}, {"Ca", 40.078}, {"Mn", 54.938045},
    {"Fe", 55.845}, {"Cu", 63.546}, {"Zn", 65.38}, {"Ga", 69.723},
    {"Ge", 72.63}, {"As", 74.9216}, {"Se", 78.96}, {"Br", 79.904},
    {"Ag", 107.8682} /* silver */, {"Au", 196.966569} /* gold */,
    {"Hg", 200.59} /* mercury */,
    {NULL, 0}};
  int k;

  for (k = 0; elemtab[k].name != NULL; k++)
    if (strcmp(elem, elemtab[k].name) == 0)
      return (real) elemtab[k].mass;
  fprintf(stderr, "unknown element [%s]\n", elem);
  return 12.0f;
}

static rv3_t *load(const char *fn, int *pn, pdbmodel_t **pmdl, real **pw)
{
  pdbmodel_t *m;
  pdbaac_t *c = NULL;
  int i, n, nres = 0;
  real (*x)[3], *w;

  die_if ((m = pdbm_read(fn, verbose)) == NULL,
      "cannot read pdb %s\n", fn);
  if (aligntype >= ALN_CA && aligntype <= ALN_MC) {
    die_if ((c = pdbaac_parse(m, verbose)) == NULL,
      "failed to parse the pdb %s\n", fn);
    nres = c->nres;
  }
  if (aligntype == ALN_CA) {
    n = nres;
    xnew(x, n);
    for (i = 0; i < n; i++)
     rv3_copy(x[i],   c->res[i].xca);
    w = NULL;
  } else if (aligntype == ALN_BB) {
    n = nres*3;
    xnew(x, n);
    xnew(w, n);
    for (i = 0; i < nres; i++) {
      rv3_copy(x[3*i],   c->res[i].xn);
      rv3_copy(x[3*i+1], c->res[i].xca);
      rv3_copy(x[3*i+2], c->res[i].xc);
      w[3*i] = 14;
      w[3*i+1] = 12;
      w[3*i+2] = 12;
    }
  } else if (aligntype == ALN_MC) {
    n = nres*4;
    xnew(x, n);
    xnew(w, n);
    for (i = 0; i < nres; i++) {
      rv3_copy(x[4*i],   c->res[i].xn);
      rv3_copy(x[4*i+1], c->res[i].xca);
      rv3_copy(x[4*i+2], c->res[i].xc);
      rv3_copy(x[4*i+3], c->res[i].xo);
      w[4*i] = 14;
      w[4*i+1] = 12;
      w[4*i+2] = 12;
      w[4*i+3] = 16;
    }
  } else if (aligntype == ALN_ALL || aligntype == ALN_HEAVY) {
    int nhvy = 0;
    x = m->x; /* use the buffer in pm */
    n = m->natm;
    xnew(w, n);
    for (i = 0; i < n; i++) {
      w[i] = elem2mass(m->atm[i].elem);
      if (aligntype == ALN_HEAVY &&
         (strcmp(m->atm[i].elem, "H") == 0 || strcmp(m->atm[i].elem, "D") == 0))
        w[i] = 0.f; /* annihilate hydrogen atoms */
      else nhvy++;
    }
    printf("%d heavy atoms / %d all atoms\n", nhvy, n);
  } else fatal("unknown alignment type %d\n", aligntype);
  *pmdl = m;
  *pn = n;
  *pw = w;
  if (c) pdbaac_free(c);
  return x;
}

static void pdbm_rottrans(pdbmodel_t *m, real rot[3][3], real trans[3])
{
  real x1[3];
  int i;

  for (i = 0; i < m->natm; i++) {
    rm3_mulvec(x1, rot, m->atm[i].x);
    rv3_inc(x1, trans);
    rv3_copy(m->atm[i].x, x1);
  }
}

int main(int argc, char **argv)
{
  pdbmodel_t *m1, *m2;
  real *w1, *w2, rmsd;
  rv3_t *x1, *x2;
  real rot[3][3], trans[3];
  int n1, n2, n;

  doargs(argc, argv);

  x1 = load(fn1, &n1, &m1, &w1);
  x2 = load(fn2, &n2, &m2, &w2);
  n = (n1 < n2) ? n1 : n2;
  if (n1 != n2) {
    fprintf(stderr, "Warning: # of atoms mismatch %s %d != %s %d, use %d\n", fn1, n1, fn2, n2, n);
  }
  rmsd = rv3_rmsd(x2, NULL, x1, w2, n, rot, trans);
  printf("rmsd = %g A\n", rmsd);
  if (verbose) {
    rm3_print(rot, "Rotation", "%8.3f", 1);
    rv3_print(trans, "Translation", "%8.3f", 1);
  }
  pdbm_rottrans(m2, rot, trans);
  printf("new %s --> %s\n", fn2, fnout);
  pdbm_write(m2, fnout);

  pdbm_free(m1);
  pdbm_free(m2);
  if (w1) free(w1);
  if (w2) free(w2);
  if (aligntype <= 1) {
    free(x1);  free(x2);
  }
  return 0;
}

