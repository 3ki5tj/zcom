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

#if defined(Macintosh) || defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

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

float dotsz = 0.002; /* size of a dot */
int width = 800, height = 400; /* preferred window size */
float xmargin = 0.1, ymargin = 0.05;
float aspect = 1.0; /* y/x for the main panel */
float eheight = 0.3; /* height of energy histogram */
float dhdescale = 2.0; /* vertical scale for dhde curve */
int pause;
double ehis[ECNT];

/* colors */
GLfloat cDhde[] = {1.0f, 0.0f, 0.0f, 0.f};  /* dhde curve */
GLfloat cDhdeGrid[] = {0.3f, 0.0f, 0.0f, 0.f}; /* dhde grids */

GLfloat cDark[] = {0.0f, 0.0f, 0.25f, 0.f}; /* background energy interval bands */
GLfloat cDarker[] = {0.0f, 0.0f, 0.15f, 0.f};

GLfloat cGray[] = {0.6f, 0.6f, 0.6f, 0.f};
GLfloat cCurr[] = {0.1f, 1.0f, 0.2f, 0.f};

GLfloat cEhis1[] = {0.5f, 0.5f, 0.0f, 0.f}; /* energy histogram bright color */
GLfloat cEhis2[] = {0.3f, 0.1f, 0.0f, 0.f}; /* energy histogram dark color */


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
    ehis[potts->E - EMIN] += 1.0;
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

/* reshape function */
static void reshape(int w, int h)
{
  float asp = 1.0f * h / w; /* screen aspect ratio */
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  aspect = asp;
  glOrtho(-xmargin, 1.f + xmargin, 
    -ymargin*asp, (1 + eheight + ymargin)*asp, -5.0f, 5.0f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/* set the current drawing color */
static void setColor(GLfloat c[]) { glColor3fv(c); }

/* convert dhde value to vertical screen scale */
static double dhde2y(double x) { return ((x - 1.0) * dhdescale + 0.5) * aspect; }

/* draw energy histogram */
static void drawEHis(void)
{
  int i, ie;
  double hmax = 0., x, y, y0 = aspect;

  for (i = 0; i < ECNT; i++)
    if (ehis[i] > hmax) hmax = ehis[i];
  setColor(cEhis2);
  for (ie = EMIN; ie < EMAX; ie += EDEL) {
    glBegin(GL_POLYGON);
    x = (ie - tmh->erg0)/(tmh->erg1 - tmh->erg0);
    y = eheight * ehis[ie - EMIN] / hmax;
    glVertex2d(x, y0);
    setColor(cEhis1);
    glVertex2d(x, y0 + aspect * y);
    x = (ie + EDEL - tmh->erg0)/(tmh->erg1 - tmh->erg0);
    glVertex2d(x, y0 + aspect * y);
    setColor(cEhis2);
    glVertex2d(x, y0);
    glEnd();
  }
}

/* draw a square dot (for a footprint) */
static void drawDot(double x, double y)
{
  double dotszy = dotsz;
  y *= aspect;
  glBegin(GL_POLYGON);
  glVertex2d(x - dotsz, y - dotszy);
  glVertex2d(x - dotsz, y + dotszy);
  glVertex2d(x + dotsz, y + dotszy);
  glVertex2d(x + dotsz, y - dotszy);
  glEnd();
}

/* draw the current position and some past point */
static void drawFootprint(double x, double y, double w, double h)
{
#define XYCNT 1024
  static double xy[XYCNT][2]; /* history of trace */
  static int xyfull = 0, xyid = 0;
  GLfloat c[3];
  int j, i, ic, icnt;

  /* add the current point to the queue */
  xy[xyid][0] = x;
  xy[xyid][1] = y;
  xyid = (xyid + 1) % XYCNT;
  if (xyid == 0) xyfull = 1;

  /* draw past points */
  icnt = xyfull ? XYCNT : xyid;
  for (ic = icnt - 1; ic >= 0; ic--) {
    i = (xyid - ic + XYCNT) % XYCNT; /* actual index */
    for (j = 0; j < 3; j++) /* past points dimmer */
      c[j] = cCurr[j]*(1.f - 1.f*ic/icnt);
    setColor(c);
    drawDot(xy[i][0], xy[i][1]);
  }

  /* draw the current point */
  y *= aspect;
  setColor(cCurr);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  glVertex2d(x - w, y); glVertex2d(x + w, y);
  glVertex2d(x, y - h); glVertex2d(x, y + h);
  glVertex2d(x - w*.3, y - h); glVertex2d(x + w*.3, y - h);
  glVertex2d(x - w*.3, y + h); glVertex2d(x + w*.3, y + h);
  glVertex2d(x - w, y - h*.3); glVertex2d(x - w, y + h*.3);
  glVertex2d(x + w, y - h*.3); glVertex2d(x + w, y + h*.3);
  glEnd();
}

/* display */
static void display(void)
{
  int i;
  double x, y;

  if (pause) return;

  /* run a few steps of tmh */
  tmhrun(speed);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  /* energy background strips */
  for (i = -1; i <= tmh->ergn; i++) {
    setColor(i % 2 == 0 ? cDark : cDarker);
    glBegin(GL_POLYGON);
    x = (i >= 0) ? (1.0*i/tmh->ergn) : -1.0;
    glVertex2d(x, 0);
    glVertex2d(x, aspect);
    x = (i < tmh->ergn) ? (1.0*(i+1)/tmh->ergn) : 2.0;
    glVertex2d(x, aspect);
    glVertex2d(x, 0);
    glEnd();
  }

  /* draw energy histogram */
  drawEHis();

  /* axes */
  setColor(cGray);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  glVertex2d(0, 0); glVertex2d(1, 0); /* x-axis */
  glVertex2d(0, aspect); glVertex2d(1, aspect);
  glVertex2d(0, 0); glVertex2d(0, aspect); /* y-axis */
  glVertex2d(1, 0); glVertex2d(1, aspect);
  glEnd();

  /* reference E-T */
  setColor(cCurr);
  glLineWidth(1.0f);
  glBegin(GL_LINES);
  glVertex2d(0, 0); glVertex2d(1.0, aspect);
  glEnd();

  /* dhde grid */
  setColor(cDhdeGrid);
  glLineStipple(1, 0x5555);
  glEnable(GL_LINE_STIPPLE);
  glBegin(GL_LINES);
  for (x = 0; x < 10.0; x += 0.1) {
    y = dhde2y(x);
    if (y < 0 || y > aspect) continue;
    glVertex2d(0, y);
    glVertex2d(1.0, y);
  }
  glEnd();
  glDisable(GL_LINE_STIPPLE);

  /* dhde reference 1.0 */
  setColor(cDhde);
  glLineWidth(1.0f);
  glBegin(GL_LINES);
  glVertex2d(0, 0.5 * aspect); glVertex2d(1.0, .5 * aspect);
  glEnd();

  /* dhde */
  glLineWidth(2.0);
  glBegin(GL_LINE_STRIP);
  for (i = 0; i <= tmh->ergn; i++) {
    x = 1.0*i/tmh->ergn;
    y = dhde2y(tmh->dhde[i]);
    glVertex2d(x, y);
  }
  glEnd();
  glLineWidth(1.0);

  /* draw trace */
  x = (potts->E - tmh->erg0)/(tmh->erg1 - tmh->erg0);
  y = (tmh->tp - tmh->tp0)/(tmh->tp1 - tmh->tp0);
  drawFootprint(x, y, 0.02*aspect, 0.03*aspect);

  glutSwapBuffers();
}


/* regularly redisplay */
static void timer(int value)
{
  glutPostRedisplay();
  glutTimerFunc(delay, timer, ++value);
}

static void keyboard(unsigned char key, int x, int y)
{
  (void) x; (void) y;
  if (key == 27 || key == 'q') {
    tmhexit(); 
  } else if (key == 'p' || key == ' ') {
    pause = !pause;
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
  tmh_initwlcvg(tmh, ampc, ampmax, sqrt(0.1), 0.95, 0, 0, 0);
  printf("erange (%g, %g), active (%g, %g)\n", tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);

  tpinit = .5 * (tp0 + tp1);
  mcrun(tmh, potts, 1.0/tpinit, tmcrun); /* equilibration */
  tmh_settp(tmh, tpinit);

  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow("TMH demo");
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutDisplayFunc(display);
  glutTimerFunc(delay, timer, 0);
  glutMainLoop();

  return 0;
}
