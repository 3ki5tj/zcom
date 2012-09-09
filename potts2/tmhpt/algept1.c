/* Algrihtm E on Potts model, zeroth/first order dhde, without using TMH */
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
#include "alge.h"

int ergmin = -1760, ergmax = -832, derg = 16, dhdeorder = 0;
double trun = 1000000*1e4;
double amp0 = 1e-2, ampc = 30.0, beta0 = 1.4;

/* configuration move under modified Hamiltonian, return de */
static int move(potts_t *pt, algei_t *al)
{
  static int steps = 0;
  int id, so, sn, eo, en, de, nb[PT2_Q], acc, ie = -1, out = 0;
  double dh, alf;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  ie = (pt->E - al->emin)/al->edel;
  die_if (ie < 0 || ie >= al->n, "ie %d out of range (0, %d)\n", ie, al->n);
  de = nb[so] - nb[sn];
  if (de != 0) {
    eo = pt->E;
    en = eo + de;
    dh = al->order ? algei_dh1(al, en, eo) : algei_dh0(al, en, eo);
    //dh = de*al->dh[ie]; /* quick and dirty way for dh  */
    out = (en >= al->emax || en < al->emin);
    acc = (dh <= 0 || rnd0() < exp(-dh));
  } else acc = 1;

  //al->cnt[ie] += 1;
  if (acc && de != 0) { /* update dh */
    if (al->order > 0)
      ie = (int)( (pt->E - al->emin + 0.499999*al->edel)/al->edel );
    al->cnt[ie] += 1;
    alf = dblmin(al->a0, al->ac/al->cnt[ie]);
    if (alf < 1e-4) alf = 1e-4;
    //if (fmod(al->cnt[ie], 1e5) < 1.0) printf("ie %d, cnt %g, amp %g\n", ie, (double) al->cnt[ie], alf);      
    al->dh[ie] += alf*de;
    /* old code for first order */
    //static int step;
    //al->dh[ie]   += alf*de*.5;
    //al->dh[ie+1] += alf*de*.5;
  }
  if (acc && !out) {
    PT2_FLIP(pt, id, so, sn, nb);
  }
  return acc * de;
}

/* change a random site */
static void randomflip(potts_t *pt)
{
  int id, so, sn, nb[PT2_Q];
  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  PT2_FLIP(pt, id, so, sn, nb);
}

static int savehis(double *his, double *acc, double *inc, const char *fn)
{
  int i;
  FILE *fp;

  xfopen(fp, fn, "w", return -1);
  for (i = 0; i < ECNT; i++) {
    if (his[i] < 0.5)  continue;
    fprintf(fp, "%d %g %g %g\n", EMIN + i * EDEL, his[i], acc[i], inc[i]);
  }
  fclose(fp);
  return 0;
}

/* entropic sampling */
static int run(potts_t *pt, double trun)
{
  algei_t *al;
  double *ehis, *eacc, *einc, it;
  int ie, de;
  
  al = algei_open(ergmin, ergmax, derg, amp0, ampc, beta0, dhdeorder);
  /* equilibrate the system till pt->E > ergmin */
  for (; pt->E <= ergmin + 4;) {
    randomflip(pt);
  }
  printf("equilibration finished, E %d\n", pt->E);
  
  xnew(ehis, ECNT);
  xnew(eacc, ECNT);
  xnew(einc, ECNT);
  for (it = 1; it <= trun; it++) {
    ie = (pt->E - EMIN)/EDEL;
    ehis[ie] += 1;
    de = move(pt, al);
    eacc[ie] += (de != 0);
    einc[ie] += abs(de);
  }
  savehis(ehis, eacc, einc, "algei.ehis");
  algei_save(al, "algei.e");
  algei_close(al);
  free(ehis); free(eacc); free(einc);
  return 0;
}

int main(void)
{
  potts_t *pt = pt2_open(L, PT2_Q);
 
  run(pt, trun);
  pt2_close(pt);
  mtsave(NULL);
  return 0;
}
