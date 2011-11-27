/* graphics version: simplified tmhab, only 3D */
#define HAVE_REAL 1
/* typedef float real; */
typedef double real;

#include "time.h"
#define ZCOM_PICK
#define ZCOM_CFG
#define ZCOM_ARGOPT
#define ZCOM_LOG
#define ZCOM_ABPRO
#define ZCOM_AV
#define ZCOM_TMH
#define ZCOM_GLEZ
#include "zcom.h"
#include "gtmh.h"

const char *fntp = "tmhab.t", *fndhde = "tmhab.e", *fnehis = "tmhab.ehis";
const char *fnpos = "ab.pos";
const char *fncfg = "tmhab.cfg";
const char *fnlog = "tmhab.tr";
double nsteps = 1000000*1000;
int seqid = 10, d = 3, model = 2;
int tmh_dhdeorder = 1;
double tmh_dhdemin = 0.1, tmh_dhdemax = 10.0;
real mddt = 5e-3f;
real thermdt = 5e-3f;
int lgclevel = 0; /* local geometry optimization: aba */
unsigned milcshake = AB_MILCSHAKE;

double tmh_ensexp = 0.0;
double tmh_emin, tmh_emax, tmh_de = 1.0;
double tmh_tps = 0.3; /* thermostat temperature */
double tmh_tp0 = 0.1, tmh_tp1 = 1.0, tmh_dtp = 0.0;
double tmh_erg0, tmh_erg1, tmh_derg = 1.0;
double tmh_elimit = 1e9;
double tmh_springk = 0.0;
double tmh_ampmax = 1e-4, tmh_ampc = 2.0;
double tmh_entampmax = 1e-2, tmh_entampc = 100.0;
int tmh_updampc = 0;
double tmh_lgvdt = 2e-3;

/* parameters for guess the erange */
int tmh_teql = 2000;
int tmh_entropic = 0; /* entropic sampling */

double nsttrace = 1e4;
double nstsave = 1e6;
int nstmv = 10;

double maxtime = 1e9;

int tmh_srand = 0; /* seed for initial random configuration */

int verbose = 1;

#define CFGGETI(x) { cfg_add(cfg, #x, "%d",  &x, #x); }
#define CFGGETD(x) { cfg_add(cfg, #x, "%lf", &x, #x); }
#define CFGGETR(x) { cfg_add(cfg, #x, "%r",  &x, #x); }

/* load setting from configuration */
static int loadcfg(const char *fn)
{
  cfg_t *cfg;

  die_if ((cfg = cfgopen(fn)) == NULL, "cannot open %s\n", fn);

  cfg_add(cfg, "seqid", "!%d", &seqid, "sequence id");
  cfg_add(cfg, "model", "!%d", &model, "model (1 or 2)");

  CFGGETD(nsteps);
  CFGGETR(mddt);
  CFGGETR(thermdt);
  CFGGETI(lgclevel);

  CFGGETD(tmh_tp0);
  CFGGETD(tmh_tp1);
  CFGGETD(tmh_dtp);

  CFGGETD(tmh_de);
  CFGGETD(tmh_derg);

  CFGGETD(tmh_tps);
  
  cfg_add(cfg, "tmh_emin", "!%lf", &tmh_emin, "energy min");
  cfg_add(cfg, "tmh_emax", "!%lf", &tmh_emax, "energy max");
  cfg_add(cfg, "tmh_erg0", "!%lf", &tmh_erg0, "active energy 0");
  cfg_add(cfg, "tmh_erg1", "!%lf", &tmh_erg1, "active energy 1");
  CFGGETI(tmh_teql);

  CFGGETD(tmh_elimit);
  CFGGETD(tmh_springk);

  CFGGETI(tmh_dhdeorder);
  CFGGETD(tmh_dhdemin);
  CFGGETD(tmh_dhdemax);
  CFGGETD(tmh_ensexp);
  CFGGETD(tmh_ampmax);
  CFGGETD(tmh_ampc);
  CFGGETI(tmh_updampc);
  CFGGETD(tmh_lgvdt);
  
  CFGGETI(tmh_entropic);
  CFGGETD(tmh_entampc);
  CFGGETD(tmh_entampmax);

  CFGGETI(tmh_srand);

  CFGGETD(nsttrace);
  CFGGETD(nstsave);
  cfg_add(cfg, "nstmv", "%d", &nstmv, "number of steps per tempering");

  die_if (cfg_match(cfg, CFG_VERBOSE|CFG_CHECKUSE) != 0, "failed match %d\n", fn);
  if (tmh_srand) srand((unsigned) tmh_srand);
  cfg_close(cfg);
  return 0;
}

/* run constant temperature run */
static void ctrun(abpro_t *ab, double tp, int teql)
{
  int t;

  for (t = 1; t <= teql; t++) {
    ab_vv(ab, (real)( tmh_tps/tp ), (real)mddt, AB_SOFTFORCE|milcshake);
    ab_rmcom(ab, ab->x, ab->v);
    ab_vrescale(ab, (real) tmh_tps, (real)thermdt);
    if (t % 100 == 0 && ab->lgcon && ab->lgact < ab->lgcnt)
      ab_updconstr(ab, 0);
  }
}

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(ARGOPT_LONGOPT); /* for -h */
  argopt_regopt(ao, "-cfg", NULL, &fncfg, "configuration file");
  argopt_regopt(ao, "-g", NULL, &fnlog, "log file");
  argopt_reghelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

abpro_t *ab;
tmh_t *tmh;
logfile_t *lg;

int delay = 100;
int speed = 200;
int pause = 0;
double t = 0.0;
int mainwinid;
real ballfudge = 0.8f;

int tmhw = 600, tmhh = 300;
float dotsz = 0.002; /* size of a dot */
float xmargin = 0.1, ymargin = 0.05;
float aspect = 1.0; /* y/x for the main panel */
float eheight = 0.3; /* height of energy histogram */
float dhdescale = 2.0; /* vertical scale for dhde curve */

/* last number is the radius */
GLfloat stcolorrad[5] = {0.8f, 0.8f, 0.8f, 1.f, .04};
GLfloat aacolorrad[2][5] = {
  {1.0f, 0.0f, 0.0f, 1.0f, 0.3f}, /* A */
  {0.0f, 0.0f, 1.0f, 1.0f, 0.3f}, /* B */
};

/* tmh run */
static void timer(int ival)
{
  int it = 0, i;
  double dhde;

  for (i = 0; i < speed; i++, t++) {
    dhde = tmh_getdhde(tmh, ab->epot) * tmh_tps / tmh->tp;
    ab_vv(ab, (real) dhde, (real) mddt, AB_SOFTFORCE|milcshake);
    ab_rmcom(ab, ab->x, ab->v);
    ab_vrescale(ab, (real) tmh_tps, (real) thermdt);
    
    if (ab->epot < ab->emin + 0.05 || (tmh->iec < 3 && rnd0() < 1e-4)) {
      double em = ab->emin;
      int lgcon = ab->lgcon;
      ab->lgcon = 0; /* turn off LGC during energy minimization */
      if (ab_localmin(ab, ab->x, 0, 0., 0, 0., AB_LMREGISTER|AB_LMWRITE) < em)
        printf("emin = %10.6f from %10.6f t %g tp %g.%30s\n", ab->emin, ab->epot, t, tmh->tp, "");
      ab->lgcon = lgcon;
    }
    gtmh_eadd(ab->epot);

    if (++it % nstmv == 0) {
      it = 0;
      if (ab->lgcon && ab->lgact < ab->lgcnt)
        ab_updconstr(ab, 0);
      if (tmh_entropic) {
        tmh_ezmoves(tmh, ab->epot, 1.0);
      } else { /* tweak updating factor by tps/tp */
        tmh_ezmove(tmh, ab->epot, tmh_tps/tmh->tp, tmh_lgvdt);
      }
    }
  
    if ((int) fmod(t, nsttrace) == 0) {
      log_printf(lg, "%g %g %g %d %d %g\n",
          t, ab->epot, tmh->tp, tmh->iec, tmh->itp, dhde);
      if (verbose) fprintf(stderr, "t = %g epot = %g, tp = %g, iec %d, itp %d, dhde %g;%20s\r",
          t, ab->epot, tmh->tp, tmh->iec, tmh->itp, dhde, "");
    }
    if ((int) fmod(t, nstsave) == 0) {
      mtsave(NULL);
      tmh_save(tmh, fntp, fnehis, fndhde, tmh->wl->lnf, t);
      ab_writepos(ab, ab->x, ab->v, fnpos);
    }
  }

  glutSetWindow(mainwinid);
  glutPostRedisplay();
  glutSetWindow(gtmh_winid);
  glutPostRedisplay();
  glutTimerFunc(delay, timer, ++ival);  
}

/* display function */
static void display(void)
{
  int i, j, aa, n = ab->n, d = ab->d;
  real xi[3], xj[3], xm = 0, sm, rb, rs; 

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  for (i = 0; i < n; i++) 
    for (j = 0; j < d; j++)
      xm += fabs(ab->x[i*d + j]);
  xm *= 2.0/(3 * ab->n); /* average box size */
  sm = 0.5f/xm;
  rs = stcolorrad[4]*sm;
  
  /* draw sticks */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, stcolorrad);
  for (i = 1; i < n; i++) {
    rv3_smul(rv3_copy(xj, ab->x + (i-1)*d), sm);
    rv3_smul(rv3_copy(xi, ab->x + i*d), sm);
    glez_drawstick(xj, xi, rs, 12);
  }

  /* draw balls */
  for (i = 0; i < n; i++) {
    aa = ab->type[i];
    rb = aacolorrad[aa][4]*sm*ballfudge;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, aacolorrad[aa]);
    glPushMatrix();
    rv3_smul(rv3_copy(xi, ab->x + i*d), sm);
    glTranslated(xi[0], xi[1], xi[2]);
    glutSolidSphere(rb, 16, 16);
    glPopMatrix();
  }
  glutSwapBuffers();
}

/* save everything and quit */
static void quit(void)
{
  log_close(lg);
  mtsave(NULL);
  tmh_save(tmh, fntp, fnehis, fndhde, tmh->wl->lnf, t);
  tmh_close(tmh);
  ab_writepos(ab, ab->x, ab->v, fnpos);
  ab_close(ab);
  exit(0);
}

enum {SPEED_MINUS, SPEED_PLUS,
  BALL_MINUS, BALL_PLUS, STICK_MINUS, STICK_PLUS, QUIT};

glez_menukey_t menukey_ctrl[] = {
  {SPEED_MINUS, 'p', "- speed"},
  {SPEED_PLUS,  'P', "+ speed"},
  {-1, '\0', NULL}};

glez_menukey_t menukey_look[] = {
  {BALL_MINUS,  'b', "- ball size"},
  {BALL_PLUS,   'B', "+ ball size"},
  {STICK_MINUS, 'k', "- stick size"},
  {STICK_PLUS,  'K', "+ stick size"},
  {-1, '\0', NULL}};

glez_menukey_t menukey[] = {
  {0, 0, "Control", menukey_ctrl},
  {0, 0, "Look",    menukey_look},
  {QUIT, 'q', "Quit"},
  {-1, '\0', NULL}};

/* handle menu */
static void menu(int id)
{
  real s1;

  if (id == SPEED_MINUS || id == SPEED_PLUS) {
    int sp1 = speed + (id == SPEED_PLUS ? 10 : -10);
    if (sp1 >= 10 && sp1 <= 1000000) speed = sp1;
  } else if (id == BALL_MINUS || id == BALL_PLUS) {
    s1 = ballfudge + (id == BALL_PLUS ? 0.02 : -0.02);
    if (s1 >= 0.05) ballfudge = s1;
  } else if (id == STICK_MINUS || id == STICK_PLUS) {
    s1 = stcolorrad[4] + (id == STICK_PLUS ? 0.005 : -0.005);
    if (s1 >= 0.02) stcolorrad[4] = (GLfloat) s1;
  } else if (id == QUIT) {
    quit();
  }
  printf("t %g, speed %d, ball %g, stick %g\n", t, speed, ballfudge, stcolorrad[4]);
}

/* inititial lights */
static void initlights(void)
{
  GLfloat light_position0[] = {1.0f, 2.0f, 3.0f, 0.0f};  /* Infinite light location. */
  GLfloat bg[] = {0.5f, 0.5f, 0.5f, 1.0f};

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, bg);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
}

int main(int argc, char **argv)
{
  double step0 = 0, amp, dtp;
  int isctn; 

  doargs(argc, argv);
  die_if (loadcfg(fncfg) != 0, "cannot load %s\n", fncfg);
  isctn = fexists(fndhde); /* continue if fndhde */
  ab = ab_open(seqid, d, model, 0.1);
  if (lgclevel) {
    ab_initconstr(ab, lgclevel);
    milcshake = 0;
  }
  if (isctn) { /* load previous data */
    if (0 != ab_readpos(ab, ab->x, ab->v, fnpos))
      return -1;
    if (0 != tmh_loaderange(fndhde, &tmh_tp0, &tmh_tp1, &dtp, 
          &tmh_erg0, &tmh_erg1, &tmh_derg, &tmh_emin, &tmh_emax, &tmh_de,
          &tmh_ensexp, &tmh_dhdeorder))
      return -1;
  } else {
    ctrun(ab, tmh_tps, tmh_teql);
    ab->epot = ab_energy(ab, ab->x, 0);
    printf("finish equilibration, epot %g, ekin %g\n", ab->epot, ab->ekin);
  }

  tmh = tmh_open(tmh_tp0, tmh_tp1, tmh_dtp, tmh_erg0, tmh_erg1, tmh_derg, 
      tmh_emin, tmh_emax, tmh_de, tmh_ensexp, tmh_dhdeorder);
  tmh->elimit = tmh_elimit;
  tmh->springk = tmh_springk;
  tmh_initwlcvg(tmh, tmh_entropic ? tmh_entampc : tmh_ampc, 
      tmh_entropic ? tmh_entampmax : tmh_ampmax, sqrt(0.1), 0.0, 
      0.0, tmh_updampc ? WLCVG_UPDLNFC : 0);
  printf("erange (%g, %g), active (%g, %g)\n", tmh->emin, tmh->emax, tmh->erg0, tmh->erg1);

  tmh_settp(tmh, tmh_tps); /* tp = tps */
  die_if (isctn && tmh_load(tmh, fnehis, fndhde, &amp, &step0) != 0,
    "cannot load tmh, %s, %s\n", fndhde, fnehis);
  t = step0;

  lg = log_open(fnlog);
  if (isctn) lg->flags = LOG_APPEND;

  /* create the main window */
  mainwinid = glezInitWindow(&argc, argv, 600, 600, "AB model");
  glezMenuKeyFunc(menu, NULL, menukey);
  initlights();
  glutDisplayFunc(display);

  /* create a tmh window */
  gtmh_init(tmh, tmhw, tmhh, "TMH");

  glutTimerFunc(delay, timer, 0);
  glutMainLoop();

  return 0;
}

