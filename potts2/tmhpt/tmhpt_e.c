/* algorihtm E on Potts model, zeroth and first order */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L)
#define EMAX 0
#define EDEL 1

#define ZCOM_PICK
#define ZCOM_POTTS2
#define ZCOM_TMH
#define TMH_NOCHECK  /* dangerous macro */
#include "zcom.h"

int dhdeorder = 0;
double erg0 = -1760, erg1 = -832, derg = 16;
double trun = 1000000*1000, trep = 100000;
double entampmax = 1e-2, entampc = 10.0;
const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis";

/* energy move under modified Hamiltonian */
static double move(tmh_t *m, potts_t *pt)
{
  int id, so, sn, eo, en, de, nb[PT2_Q], acc, iec, out = 0;
  double dh, amp;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  iec = (int)((pt->E - m->erg0)/m->derg);
  if (de != 0) { /* avoid calling hdif() if possible */
    eo = pt->E;
    en = eo + de;
    //dh = m->dhde[ iec ] * de;
    dh = tmh_hdif(m, en, eo); /* modified Hamiltonian */
    out = (en >= m->erg1 || en < m->erg0);
    acc = (dh <= 0) || (rnd0() < exp(-dh));
  } else acc = 1;

  // m->dhdecnt[iec] += 1; 
  if (acc && de != 0) {
    m->dhdecnt[iec] += 1;
    amp = dblmin(entampc/m->dhdecnt[iec], entampmax);
    m->dhde[iec] = dblmax(m->dhde[iec] + de * amp, 0);
  }  

  if (acc && !out) {
    PT2_FLIP(pt, id, so, sn, nb);
  }
  m->tpehis[pt->E - EMIN] += 1;
  return amp;
}

/* entropic sampling */
static int run(tmh_t *m, potts_t *pt, double trun)
{
  double t, amp;

  for (t = 1; pt->E < erg0 + 4; t++) { /* equilibration */
    int id, nb[PT2_Q], so, sn;
    PT2_PICK(pt, id, nb);
    PT2_NEWFACE(pt, id, so, sn);
    PT2_FLIP(pt, id, so, sn, nb);
  }
  printf("after equilibration erg %d\n", pt->E);

  for (t = 1; t <= trun; t += 1) { /* production */
    amp = move(m, pt);
  }
  tmh_save(m, fntp, fnehis, fndhde, amp, t);
  return 0;
}

int main(void)
{
  potts_t *pt = pt2_open(L, PT2_Q);
  tmh_t *m = tmh_open(1.0, 0.9, -0.01, erg0, erg1, derg, EMIN, EMAX + EDEL, EDEL, 2.0, dhdeorder);
  
  run(m, pt, trun);

  pt2_close(pt);
  tmh_close(m);
  mtsave(NULL);
  return 0;
}

