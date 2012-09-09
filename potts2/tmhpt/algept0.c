/* algrihtm E on Potts model, quick and dirty algorithm, without TMH */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L)
#define EMAX 0
#define EDEL 1
#define ECNT ((EMAX - EMIN)/EDEL + 1)

#define ZCOM_PICK
#define ZCOM_POTTS2
#define ZCOM_HIST
#include "zcom.h"

int ergmin = -2000, ergmax = -800, ergdel = 4;
double beta0 = 1.4, alpha = 1e-3;
double trun = 1000000*100;

typedef struct {
  int emin, emax, edel, n;
  double *dh; /* dH/dE */
} algei_t;

static algei_t *algei_open(int emin, int emax, int edel)
{
  algei_t *al;
  int i;
  
  xnew(al, 1);
  al->emin = emin;
  al->emax = emax;
  al->edel = edel;
  al->n = (emax - emin)/edel;
  xnew(al->dh, al->n + 1);
  for (i = 0; i <= al->n; i++) al->dh[i] = beta0;
  return al;
}

static void algei_close(algei_t *al)
{
  free(al->dh);
  free(al);
}

static int algei_savedh(algei_t *al, const char *fn)
{
  FILE *fp;
  int i;
  
  xfopen(fp, fn, "w", return -1);
  for (i = 0; i <= al->n; i++)
    fprintf(fp, "%d %.6f\n", al->emin + i * al->edel, al->dh[i]);
  fclose(fp);
  return 0;
}

/* configuration move under modified Hamiltonian, return de */
static void move(potts_t *pt, algei_t *al)
{
  static int steps = 0;
  int id, so, sn, eo, en, de, nb[PT2_Q], acc, ie = -1, out = 0;
  double dh, alf;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  if (de != 0) {
    eo = pt->E;
    en = eo + de;
    ie = (eo - al->emin)/al->edel;
    dh = de*al->dh[ie]; /* approximate dh, quick and dirty */
    out = (en >= al->emax || en < al->emin);
    acc = (dh <= 0 || rnd0() < exp(-dh));
  } else acc = 1;

  if (acc && de != 0) al->dh[ie] += alpha * de;
  if (acc && !out) {  PT2_FLIP(pt, id, so, sn, nb);  }
}

/* change a random site */
static void randflip(potts_t *pt)
{
  int id, so, sn, nb[PT2_Q];
  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  PT2_FLIP(pt, id, so, sn, nb);
}

/* entropic sampling */
static int run(potts_t *pt, double trun)
{
  algei_t *al;
  double *ehis;
  int it;
  
  al = algei_open(ergmin, ergmax, ergdel);
  /* equilibrate the system till pt->E > ergmin */
  for (it = 1; pt->E <= ergmin + 4; it++) randflip(pt);
  printf("equilibration finished in %d steps\n", it);
  
  xnew(ehis, ECNT);
  for (it = 1; it <= trun; it++) {
    move(pt, al);
    ehis[pt->E - EMIN] += 1;
  }
  histsave(ehis, 1, ECNT, EMIN, EDEL, 0, "algei.ehis");
  algei_savedh(al, "algei.e");
  algei_close(al);
  free(ehis);
  return 0;
}

int main(void)
{
  potts_t *pt = pt2_open(L, PT2_Q);
  run(pt, trun);
  pt2_close(pt);
  return 0;
}
