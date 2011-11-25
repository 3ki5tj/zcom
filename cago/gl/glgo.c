#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_AV
#define ZCOM_CAGO
#define ZCOM_GLEZ
#include "zcom.h"

const char *fnpdb = "pdb/1MBN.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 6.f;
cago_t *go;

real tp = 0.3f;
int verbose = 0;
int step = 0;
real mddt = 0.002f, thermdt = 0.01f;
av_t avep[1], avtk[1], avrmsd[1];

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_regarg(ao, NULL, &fnpdb, "pdbfile");
  argopt_reghelp(ao, "-h");
  argopt_regopt(ao, "-T", "%r", &tp, "Temperature");
  argopt_regopt(ao, "-v", "%d", &verbose, "verbose");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

/* initialize MD */
static int initmd(void)
{
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("%s n %d, tp %.3f, epot = %g, %g (ref), rmsd = %g\n", 
      fnpdb, go->n, tp, go->epot, go->epotref, go->rmsd);
  return 0;
}

/* reset data */
static void reset(void)
{
  step = 0;
  av_clear(avep);
  av_clear(avtk);
  av_clear(avrmsd);
}

int delay = 50;
int speed = 100;
int trep = 10000;
real ballfudge = 0.8f;

/* last number is the radius */
GLfloat stcolorrad[5] = {0.8f, 0.8f, 0.8f, 1.f, .5};
GLfloat aacolorrad[20][5] = {
  {1.0f, 1.0f, 1.0f, 1.0f, 0.9f}, /* GLY */
  {0.6f, 0.6f, 0.6f, 1.0f, 1.0f}, /* ALA */
  {0.7f, 0.7f, 0.7f, 1.0f, 1.3f}, /* VAL */
  {0.8f, 0.8f, 0.8f, 1.0f, 1.4f}, /* LEU */
  {0.8f, 0.8f, 0.8f, 1.0f, 1.4f}, /* ILE */
  {0.8f, 0.8f, 0.8f, 1.0f, 1.3f}, /* PRO */
  {0.8f, 0.5f, 1.0f, 1.0f, 1.2f}, /* SER */
  {0.5f, 0.3f, 0.8f, 1.0f, 1.3f}, /* THR */
  {1.0f, 1.0f, 0.0f, 1.0f, 1.2f}, /* CYS */
  {0.6f, 0.6f, 0.0f, 1.0f, 1.4f}, /* MET */
  {0.2f, 1.0f, 0.2f, 1.0f, 1.4f}, /* ASN */
  {0.1f, 1.0f, 0.6f, 1.0f, 1.5f}, /* GLN */
  {1.0f, 0.0f, 0.0f, 1.0f, 1.4f}, /* ASP */
  {1.0f, 0.2f, 0.0f, 1.0f, 1.5f}, /* GLU */
  {0.0f, 0.4f, 1.0f, 1.0f, 1.5f}, /* LYS */
  {0.0f, 0.0f, 1.0f, 1.0f, 1.7f}, /* ARG */
  {0.0f, 0.4f, 0.8f, 1.0f, 1.6f}, /* HIS */
  {0.8f, 0.9f, 0.9f, 1.0f, 1.7f}, /* PHE */
  {1.0f, 0.9f, 0.9f, 1.0f, 1.8f}, /* TYR */
  {0.9f, 1.0f, 1.0f, 1.0f, 2.0f}, /* TRP */
};

static void timer(int ival)
{
  int i;

  for (i = 0; i < speed; i++) {
    cago_vv(go, 1.f, mddt);
    cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, tp, thermdt);
    cago_rotfit(go, go->x, NULL);
    step++;
    av_add(avep, go->epot);
    av_add(avtk, go->tkin);
    av_add(avrmsd, go->rmsd);
  } 
  glutPostRedisplay();
  glutTimerFunc(delay, timer, 0);
  
  if (step % trep == 0 && ival == 0)
    printf("t = %.3f, U/N %10.3f(%10.3f), T %6.3f(%6.3f) rmsd %8.3f(%8.3f)\n",
      step*mddt, go->epot, av_getave(avep), go->tkin, av_getave(avtk),
      go->rmsd, av_getave(avrmsd));
}

static void display(void)
{
  int i, j, aa;
  real xi[3], xj[3], xm = 0, sm, rb, rs; 

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  for (i = 0; i < go->n; i++) 
    for (j = 0; j < 3; j++)
      xm += fabs(go->x[i][j]);
  xm *= 2.0/(3*go->n); /* average box size */
  sm = 0.5f/xm;
  rs = stcolorrad[4]*sm;
  
  /* draw sticks */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, stcolorrad);
  for (i = 1; i < go->n; i++) {
    rv3_smul(rv3_copy(xj, go->x[i-1]), sm);
    rv3_smul(rv3_copy(xi, go->x[i]), sm);
    glez_drawstick(xj, xi, rs, 12);
  }

  /* draw balls */
  for (i = 0; i < go->n; i++) {
    aa = go->aa[i];
    rb = aacolorrad[aa][4]*sm*ballfudge;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, aacolorrad[aa]);
    glPushMatrix();
    rv3_smul(rv3_copy(xi, go->x[i]), sm);
    glTranslated(xi[0], xi[1], xi[2]);
    glutSolidSphere(rb, 16, 16);
    glPopMatrix();
  }
  glutSwapBuffers();
}

enum {T_MINUS, T_PLUS, SPEED_MINUS, SPEED_PLUS,
  BALL_MINUS, BALL_PLUS, STICK_MINUS, STICK_PLUS};

glez_menukey_t menukey_ctrl[] = {
  {T_MINUS,     't', "- temperature"},
  {T_PLUS,      'T', "+ temperature"},
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
  {-1, '\0', NULL}};

static void menu(int id)
{
  real s1;

  if (id == T_MINUS || id == T_PLUS) {
    real T1 =  tp+ (id == T_PLUS ? 0.02 : -0.02);
    if (T1 > 0.049) tp = T1;
    reset(); 
  } else if (id == SPEED_MINUS || id == SPEED_PLUS) {
    int sp1 = speed + (id == SPEED_PLUS ? 10 : -10);
    if (sp1 >= 10 && sp1 <= 10000) speed = sp1;
  } else if (id == BALL_MINUS || id == BALL_PLUS) {
    s1 = ballfudge + (id == BALL_PLUS ? 0.05 : -0.05);
    if (s1 >= 0.05) ballfudge = s1;
  } else if (id == STICK_MINUS || id == STICK_PLUS) {
    s1 = stcolorrad[4] + (id == STICK_PLUS ? 0.02 : -0.02);
    if (s1 >= 0.02) stcolorrad[4] = (GLfloat) s1;
  }
  printf("T %.2f, speed %d, ball %g, stick %g\n", 
      tp, speed, ballfudge, stcolorrad[4]);
}

static void initgui(void)
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
  glezInitWindow(&argc, argv, 600, 600, "C-alpha Go model");
  glezMenuKeyFunc(menu, NULL, menukey);
  glutDisplayFunc(display);
  doargs(argc, argv);
  if (initmd() != 0) exit(1);
  initgui();
  glutTimerFunc(delay, timer, 0);
  glutMainLoop();  
  return 0;
}

