#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_AV
#define ZCOM_CAGO
#include "zcom.h"
#include <GL/glut.h>

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
  argopt_regopt_help(ao, "-h");
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
real zoomscale = 0.5f;
real ballfudge = 0.8f;
enum {MS_NONE, MS_ZOOM, MS_ROTATE};
int msaction = MS_NONE, msdown, msx, msy;

/* last number is the radius */
GLfloat stcolorrad[5] = {0.8f, 0.8f, 0.8f, 1.f, .5};
GLfloat aacolorrad[20][5] = {
  {0.5f, 0.5f, 0.2f, 1.0f, 0.9f}, /* GLY */
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

/* draw a stick from a to b, with radius r and nface faces */
static void drawStick(real a[], real b[], real r, int nface)
{
  int i;
  real l, c, s, th;
  rv3_t x, y, z, p, q, u;

  rv3_diff(z, b, a);
  l = rv3_norm(z); /* stick length */
  rv3_smul(z, 1.f/l);
  
  rv3_normalize( rv3_make(x, z[1], -z[0], 0) ); /* a vector perpendicular to v */
  rv3_normalize( rv3_cross(y, z, x) );

  glBegin(GL_QUAD_STRIP);
  for (i = 0; i <= nface; i++) {
    th = 2.*M_PI*i/nface;
    c = (real) cos(th);
    s = (real) sin(th);
    rv3_lincomb2(u, x, y, c, s);
    rv3_lincomb2(p, a, u, 1.f, r);
    rv3_lincomb2(q, b, u, 1.f, r);
    glNormal3d(u[0], u[1], u[2]);
    glVertex3d(p[0], p[1], p[2]);
    glVertex3d(q[0], q[1], q[2]);
  }
  glEnd();
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
  sm = zoomscale/xm;
  rs = stcolorrad[4]*sm;
  
  /* draw sticks */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, stcolorrad);
  for (i = 1; i < go->n; i++) {
    rv3_smul(rv3_copy(xj, go->x[i-1]), sm);
    rv3_smul(rv3_copy(xi, go->x[i]), sm);
    drawStick(xj, xi, rs, 12);
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

/* menu */
enum {T_MINUS, T_PLUS, SPEED_MINUS, SPEED_PLUS,
  BALL_MINUS, BALL_PLUS, STICK_MINUS, STICK_PLUS, SCALE_MINUS, SCALE_PLUS,
  X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS, 
  FULLSCREEN, QUIT, MENULAST};

struct { int id; unsigned char key; const char *desc; } menukey[MENULAST] = {
  {T_MINUS,     't', "- temperature"},
  {T_PLUS,      'T', "+ temperature"},
  {SPEED_MINUS, 'p', "- speed"},
  {SPEED_PLUS,  'P', "+ speed"},
  {BALL_MINUS,  'b', "- ball size"},
  {BALL_PLUS,   'B', "+ ball size"},
  {STICK_MINUS, 'k', "- stick size"},
  {STICK_PLUS,  'K', "+ stick size"},
  {SCALE_MINUS, 's', "Zoom out"},
  {SCALE_PLUS,  'S', "Zoom in"},
  {X_MINUS,     'x', "Rotate around -x"},
  {X_PLUS,      'X', "Rotate around +x"},
  {Y_MINUS,     'y', "Rotate around -y"},
  {Y_PLUS,      'Y', "Rotate around +y"},
  {Z_MINUS,     'z', "Rotate around -z"},
  {Z_PLUS,      'Z', "Rotate around +z"},
  {FULLSCREEN,  'f', "Toggle full screen"},
  {QUIT,        'q', "Quit"}
};

static void menu(int id)
{
  GLfloat mat[4][4];
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
  } else if (id == SCALE_MINUS || id == SCALE_PLUS) {
    s1 = zoomscale + (id == SCALE_PLUS ? 0.02 : -0.02);
    if (s1 >= 0.02) zoomscale = s1;
    glutPostRedisplay();
  } else if (id == X_MINUS || id == X_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == X_PLUS) ? 5.f : -5.f, mat[0][0], mat[1][0], mat[2][0]);
    glutPostRedisplay();
  } else if (id == Y_MINUS || id == Y_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == Y_PLUS) ? 5.f : -5.f, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  } else if (id == Z_MINUS || id == Z_PLUS) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == Z_PLUS) ? 5.f : -5.f, mat[0][2], mat[1][2], mat[2][2]);
    glutPostRedisplay();
  } else if (id == FULLSCREEN) {
    static int full = 0, x, y, w, h;
    full = !full;
    if (full) {
      x = glutGet(GLUT_WINDOW_X);
      y = glutGet(GLUT_WINDOW_Y);
      w = glutGet(GLUT_WINDOW_WIDTH);
      h = glutGet(GLUT_WINDOW_HEIGHT);
      glutFullScreen(); 
    } else {
      glutPositionWindow(x, y);
      glutReshapeWindow(w, h);
    }
  } else if (id == QUIT) {
    exit(0);
  }
  printf("%s: T %.2f, speed %d, ball %g, stick %g\n", 
      menukey[id].desc, tp, speed, ballfudge, stcolorrad[4]);
}

static void initmenu(void)
{
  int i;
  char s[64];

  glutCreateMenu(menu);
  for (i = 0; i < MENULAST; i++) {
    sprintf(s, "%s, key: %c\n", menukey[i].desc, menukey[i].key);
    glutAddMenuEntry(s, menukey[i].id);
  }
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

static void keypress(unsigned char c, int x, int y)
{
  int i;
  
  (void) x; (void) y;
  if (c == 27) exit(0);
  for (i = 0; i < MENULAST; i++)
    if (c == menukey[i].key)
      menu(menukey[i].id);
}

static void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    msdown++;
    if (verbose >= 1) printf("mouse %d down at %d,%d; %d\n", button, msx, msy, msdown);
    if (button == 3) {
      menu(SCALE_PLUS);
    } else if (button == 4) {
      menu(SCALE_MINUS);
    }
  } else if (--msdown <= 0) {
    msdown = 0;
    msaction = MS_NONE;
  }
  msx = x;
  msy = y;
}

static void motion(int x, int y) 
{
  if (x == msx && y == msy) return;
  if (verbose >= 2) printf("mouse at %d,%d; %d\n", x, y, msdown);
  if (msdown) {
    float angx = (float)( (y - msy)*360.f/glutGet(GLUT_WINDOW_HEIGHT) );
    float angy = (float)( (x - msx)*360.f/glutGet(GLUT_WINDOW_WIDTH) );
    float mat[4][4];

    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotated(angx, mat[0][0], mat[1][0], mat[2][0]);
    glRotated(angy, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  }
  msx = x;
  msy = y;
}

static void reshape(int w, int h)
{
  double xs = 1, ys = 1;
  if (w > h) xs = 1.*w/h;
  else ys = 1.*h/w;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-xs, xs, -ys, ys, -20, 20);
  glMatrixMode(GL_MODELVIEW);
}

static void initgui(void)
{
  GLfloat light_position0[] = {1.0f, 2.0f, 3.0f, 0.0f};  /* Infinite light location. */
  GLfloat bg[] = {0.0f, 0.0f, 0.0f, 1.0f};

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, bg);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
  gluLookAt(0.0, 0.0, 5.0,  /* eye is at (0,0,5)  */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */  
}

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitWindowSize(800, 800);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow("C-alpha Go model");
  doargs(argc, argv);
  if (initmd() != 0) exit(1);
  initgui();
  initmenu();
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutKeyboardFunc(keypress);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutTimerFunc(delay, timer, 0);
  glutMainLoop();  
  return 0;
}

