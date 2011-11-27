/* Tempering with modified Hamiltonian on Potts model
 * graphics version */
#define PT2_LB 5
#define PT2_Q  10
#define L (1 << PT2_LB)
#define EMIN (-2*L*L)
#define EMAX 0
#define ECNT (EMAX - EMIN + 1)
#define EDEL 1

#define BLOCK 10 /* number of runs */

#define ZCOM_PICK
#define ZCOM_POTTS2
#define ZCOM_TMH
#include "zcom.h"
#include "gtmh.h"

potts_t *potts;
tmh_t *tmh;

const char *fntp = "tmhpt.t", *fndhde = "tmhpt.e", *fnehis = "tmhpt.ehis", *fnpos = "pt.pos";
double tp0 = 0.67, tp1 = 0.77, dtp = 0.001;  /* for bstyle == 0 */
double erg0 = -1760, erg1 = -832, derg = 4, elimit = 1e8;
double tmcrun = 1000000;
double ampmax = 1e-6, ampc = 1.0;
double lgvdt = 1e-6;
int dhdeorder = 1;
double springk = 1000.0; /* elastic force to ensure smoothness */
int easy = 1; /* use easy move */

int delay = 100;
int speed = 100000; /* number of steps per frame */
double t = 0.0;

/* energy move under modified Hamiltonian */
static int tmhmove(tmh_t *m, potts_t *pt, double beta)
{
  int id, so, sn, eo, en, de, nb[PT2_Q], acc;
  double dh;

  PT2_PICK(pt, id, nb);
  PT2_NEWFACE(pt, id, so, sn); /* so --> sn */
  de = nb[so] - nb[sn];
  if (de > 0) { /* avoid calling hdif() if possible */
    eo = pt->E;
    en = eo + de;
    dh = tmh_hdif(m, en, eo); /* modified Hamiltonian */
    acc = (dh <= 0) || (rnd0() < exp(-dh*beta));
  } else acc = 1;
  if (acc) { PT2_FLIP(pt, id, so, sn, nb); return de; }
  else return 0;
}

/* constant temperature MC run */
static void mcrun(tmh_t *m, potts_t *pt, double beta, double tmax)
{ for (; tmax >= 0; tmax--) tmhmove(m, pt, beta); }

/* run nsteps of tmh */
static void tmhrun(int nsteps)
{
  int i, it;
  double epot, del;
  tmh_t *m = tmh;

  /* let it evolve */
  for (it = 0; it < nsteps; it++, t++) {
    tmhmove(m, potts, 1.0/tmh->tp);
    gtmh_eadd(potts->E);
    if (easy) { /* use shortcut */
      tmh_ezmove(m, potts->E, 1.0, lgvdt);
    } else { /* move it manually */
      epot = potts->E;
      tmh_eadd(m, epot);
      if (fabs(epot - m->ec) < m->elimit || epot > m->erg1 || epot < m->erg0) {
        die_if (m->wl == NULL, "call tmh_initamp first, %p\n", (void *) m->wl);
        wlcvg_update(m->wl, m->tp); /* compute updating amplitude */
        del = (epot - m->ec);
#if 0
        /* WRONG: bin must be located by ec instead of epot */
        i = intconfine((int)((epot - m->erg0)/m->derg + .5), 0, m->ergn);
#endif
        i = m->iec;
        /* apply elastic force */
        if (i > 0) del += springk * (m->dhde[i-1] - m->dhde[i]);
        if (i < (m->dhdeorder ? m->ergn : m->ergn - 1))
          del += springk * (m->dhde[i+1] - m->dhde[i]);
        del *= m->wl->lnf;
        m->dhde[i] = dblconfine(m->dhde[i] + del, m->dhdemin, m->dhdemax);
      }
      tmh_lgvmove(m, epot, lgvdt);
    }
  }
}

/* exit the program */
static void tmhexit(void)
{
  tmh_save(tmh, fntp, fnehis, fndhde, tmh->wl->lnf, t);
  pt2_save(potts, fnpos);
  pt2_close(potts);
  tmh_close(tmh);
  mtsave(NULL);
  exit(0);
}

/* regularly redisplay */
static void timer(int value)
{
  /* run a few steps of tmh */
  if (!gtmh_pause) tmhrun(speed);
  glutPostRedisplay();
  glutTimerFunc(delay, timer, ++value);
}

static void keyboard(unsigned char key, int x, int y)
{
  (void) x; (void) y;
  if (key == 27 || key == 'q') {
    tmhexit(); 
  } else if (key == 'p' || key == ' ') {
    gtmh_pause = !gtmh_pause;
  } else if (key == 'h') {
    speed /= 2; if (speed < 1) speed = 1;
  } else if (key == 'd') {
    speed *= 2;
  }
  printf("speed %d, amp %g, t %g\n", speed, tmh->wl->lnf, t);
}

int main(int argc, char **argv)
{
  double tpinit, emin = EMIN, emax = EMAX, de = EDEL, ensexp = 2.0;

  potts = pt2_open(L, PT2_Q);
  tmh = tmh_open(tp0, tp1, dtp, erg0, erg1, derg, emin, emax, de, ensexp, dhdeorder);
  tmh->elimit = elimit;
  tmh->springk = springk;
  tmh_initwlcvg(tmh, ampc, ampmax, sqrt(0.1), 0.95, 0, 0);
  printf("erange (%g, %g), active (%g, %g)\n", tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);

  tpinit = .5 * (tp0 + tp1);
  mcrun(tmh, potts, 1.0/tpinit, tmcrun); /* equilibration */
  tmh_settp(tmh, tpinit);

  glutInit(&argc, argv);

  gtmh_init(tmh, 800, 400, "TMH demo");
  glutKeyboardFunc(keyboard);
  glutTimerFunc(delay, timer, 0);
  glutMainLoop();

  return 0;
}
