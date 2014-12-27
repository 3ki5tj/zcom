#include "rv3.h"
#include "util.h"
#ifndef GLEZ_H__
#define GLEZ_H__
#if defined(Macintosh) || defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



static int glez_x, glez_y; /* current position */
static int glez_msdown; /* mouse state */

typedef struct tag_glez_menukey_t {
  int id;
  int key;
  const char *str;
  struct tag_glez_menukey_t *sub; /* pointer to sub menus */
} glez_menukey_t;

static glez_menukey_t *glez_user_menukey;

enum { GLEZ_MENU0 = 1000,
  GLEZ_MOVXM, GLEZ_MOVXP, GLEZ_MOVYM, GLEZ_MOVYP, GLEZ_MOVZM, GLEZ_MOVZP,
  GLEZ_ROTXM, GLEZ_ROTXP, GLEZ_ROTYM, GLEZ_ROTYP, GLEZ_ROTZM, GLEZ_ROTZP,
  GLEZ_SCLM, GLEZ_SCLP, GLEZ_FULLS, GLEZ_MENU1 };

/* rotation sub-menu */
glez_menukey_t glez_menukey_rot[] = {
  {GLEZ_ROTXM,  'x', "Rotate around -x", NULL},
  {GLEZ_ROTXP,  'X', "Rotate around +x", NULL},
  {GLEZ_ROTYM,  'y', "Rotate around -y", NULL},
  {GLEZ_ROTYP,  'Y', "Rotate around +y", NULL},
  {GLEZ_ROTZM,  'z', "Rotate around -z", NULL},
  {GLEZ_ROTZP,  'Z', "Rotate around +z", NULL},
  {-1, '\0', NULL, NULL}};

/* scaling sub-menu */
glez_menukey_t glez_menukey_scl[] = {
  {GLEZ_SCLM,   's', "Zoom out", NULL},
  {GLEZ_SCLP,   'S', "Zoom in", NULL},
  {GLEZ_FULLS,  'f', "Toggle fullscreen", NULL},
  {-1, '\0', NULL, NULL}};

glez_menukey_t glez_menukey_mov[] = {
  {GLEZ_MOVXM,  'l', "Move toward -x", NULL},
  {GLEZ_MOVXP,  'r', "Move toward +x", NULL},
  {GLEZ_MOVYM,  'd', "Move toward -y", NULL},
  {GLEZ_MOVYP,  'u', "Move toward +y", NULL},
  {GLEZ_MOVZM,  'f', "Move toward -y", NULL},
  {GLEZ_MOVZP,  'n', "Move toward +y", NULL},
  {-1, '\0', NULL, NULL}};

/* main menu */
glez_menukey_t glez_menukey[] = {
  {0,           0, "Move",   glez_menukey_mov},
  {0,           0, "Rotate", glez_menukey_rot},
  {0,           0, "Zoom",   glez_menukey_scl},
  {-1, '\0', NULL, NULL}};

static void (*glez_user_reshapefunc)(int w, int h) = NULL;
#define glezReshapeFunc(f)  { glez_user_reshapefunc = f;  glutReshapeFunc(glez_reshapefunc); }

static void (*glez_user_menufunc)(int id) = NULL;
static void (*glez_user_keyboardfunc)(unsigned char c, int x, int y) = NULL;

static void (*glez_user_mousefunc)(int button, int state, int w, int h) = NULL;
#define glezMouseFunc(f)    { glez_user_mousefunc = f;  glutMouseFunc(glez_mousefunc); }

static void (*glez_user_motionfunc)(int w, int h) = NULL;
#define glezMotionFunc(f)   { glez_user_motionfunc = f;  glutMotionFunc(glez_motionfunc); }



#ifdef _OPENMP
#pragma omp threadprivate(glez_x, glez_y, glez_msdown)
#pragma omp threadprivate(glez_user_menukey)
#pragma omp threadprivate(glez_user_reshapefunc)
#pragma omp threadprivate(glez_user_menufunc)
#pragma omp threadprivate(glez_user_keyboardfunc)
#pragma omp threadprivate(glez_user_mousefunc)
#pragma omp threadprivate(glez_user_motionfunc)
#endif



/* toggle full screen state */
INLINE void glez_fullscreen(void)
{
  static int full = 0, x, y, w, h;
#ifdef _OPENMP
#pragma omp threadprivate(full, x, y, w, h)
#endif

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
}



/* standard reshape function for GLUT
 * x: (-1, 1), y: (-1, 1), z: (-10, 10) */
INLINE void glez_reshapefunc(int w, int h)
{
  double xs = 1, ys = 1, zs = 10;

  if (w > h) xs = 1.*w/h;
  else ys = 1.*h/w;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-xs, xs, -ys, ys, -zs, zs);
  glMatrixMode(GL_MODELVIEW);

  if (glez_user_reshapefunc) (*glez_user_reshapefunc)(w, h);
}



/* menu function */
INLINE void glez_menufunc(int id)
{
  GLfloat mat[4][4], amp;

  if (id <= GLEZ_MENU0 || id >= GLEZ_MENU1) {
    if (glez_user_menufunc) (*glez_user_menufunc)(id);
  } else if (id == GLEZ_SCLM || id == GLEZ_SCLP) { /* scaling */
    GLfloat s = (id == GLEZ_SCLP) ? 1.02f : 1.0f/1.02f;

    glScalef(s, s, s);
    glutPostRedisplay();
  } else if (id == GLEZ_FULLS) { /* full screen */
    glez_fullscreen();
    glutPostRedisplay();
  } else if (id >= GLEZ_MOVXM && id <= GLEZ_MOVZP) { /* translations */
    amp = ((id + GLEZ_MOVXM) % 2) ? .02f : -.02f;
    id = (id - GLEZ_MOVXM) / 2;
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glTranslatef(mat[0][id]*amp, mat[1][id]*amp, mat[2][id]*amp);
    glutPostRedisplay();
  } else if (id >= GLEZ_ROTXM && id <= GLEZ_ROTZP) { /* rotations */
    if (id < GLEZ_ROTXM || id > GLEZ_ROTZP) return;
    amp = ((id + GLEZ_ROTXM) % 2) ? 5.f : -5.f;
    id = (id - GLEZ_ROTXM) / 2;
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef(amp, mat[0][id], mat[1][id], mat[2][id]);
    glutPostRedisplay();
  }
}



/* recursively match short-cut keys */
INLINE int glez_keylow(unsigned char c, glez_menukey_t *mk)
{
  int i;

  if (mk == NULL) return 0;
  for (i = 0; ; i++) {
    //printf("%d %c vs %c(%d) id %d, str %s, sub %p\n", i, c, mk[i].key, mk[i].key, mk[i].id, mk[i].str, mk[i].sub);
    if (mk[i].sub != NULL) {
      if (glez_keylow(c, mk[i].sub)) return 1;
    } else if (mk[i].key == 0) {
      break;
    } else if (c == mk[i].key) {
      glez_menufunc(mk[i].id);
      return 1;
    }
  }
  return 0;
}



/* system keyboard function */
INLINE void glez_keyboardfunc(unsigned char c, int x, int y)
{
  if (glez_keylow(c, glez_menukey)) return;
  if (glez_keylow(c, glez_user_menukey)) return;
  /* pass an unhandled key to glez_user_keyboardfunc */
  if (glez_user_keyboardfunc) (*glez_user_keyboardfunc)(c, x, y);
}



/* create a menu hierarchy as specified by mk */
INLINE void glez_addmenu(glez_menukey_t *mk)
{
  int i, menuid, subid;
  char s[64];

  if (mk == NULL) return;
  for (i = 0; ; i++) {
    //printf("%d: id %d, str %s, sub %p\n", i, mk[i].id, mk[i].str, mk[i].sub);
    if (mk[i].sub != NULL) { /* sub menu */
      menuid = glutGetMenu(); /* get menu id */
      subid = glutCreateMenu(glez_menufunc);
      glez_addmenu(mk[i].sub); /* recursively create sub-menus */
      glutSetMenu(menuid); /* return to the menu */
      glutAddSubMenu(mk[i].str, subid);
    } else if (mk[i].key == 0) {
      break;
    } else { /* regular menu */
      sprintf(s, "%.32s, key: %c\n", mk[i].str, mk[i].key);
      glutAddMenuEntry(s, mk[i].id);
    }
  }
}



INLINE void glezMenuKeyFunc(void (*menuf)(int), void (*keyf)(unsigned char, int, int),
    glez_menukey_t *mk)
{
  glez_user_menufunc = menuf;
  glez_user_keyboardfunc = keyf;
  glez_user_menukey = mk;

  glutCreateMenu(glez_menufunc);
  glez_addmenu(mk); /* add user menu */
  glez_addmenu(glez_menukey); /* add system menu */
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  glutKeyboardFunc(glez_keyboardfunc);
}



INLINE void glez_mousefunc(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    glez_msdown++;
    if (button == 3)
      glez_menufunc(GLEZ_SCLP);
    else if (button == 4)
      glez_menufunc(GLEZ_SCLM);
  } else if (--glez_msdown <= 0) {
    glez_msdown = 0;
  }
  glez_x = x;
  glez_y = y;
  if (glez_user_mousefunc) (*glez_user_mousefunc)(button, state, x, y);
}



/* mouse motion function for GLUT */
INLINE void glez_motionfunc(int x, int y)
{
  if (x == glez_x && y == glez_y) return;
  if (glez_msdown) {
    float angx = (float)( (y - glez_y) * 360.f / glutGet(GLUT_WINDOW_HEIGHT) );
    float angy = (float)( (x - glez_x) * 360.f / glutGet(GLUT_WINDOW_WIDTH) );
    float mat[4][4];

    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef(angx, mat[0][0], mat[1][0], mat[2][0]);
    glRotatef(angy, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  }
  glez_x = x;
  glez_y = y;
  if (glez_user_motionfunc) (*glez_user_motionfunc)(x, y);
}



int glezInitWindow(int *argc, char **argv, int w, int h, const char *name)
{
  int winid;

  glutInit(argc, argv);
  w = intmin(w, glutGet(GLUT_SCREEN_HEIGHT));
  h = intmin(h, glutGet(GLUT_SCREEN_WIDTH));
  glutInitWindowSize(w, h);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  winid = glutCreateWindow(name);

  /* register glez functions */
  glezReshapeFunc(NULL);
  glezMenuKeyFunc(NULL, NULL, NULL);
  glezMouseFunc(NULL);
  glezMotionFunc(NULL);

  return winid;
}



/* draw a stick from a to b, with radius r and nface faces */
INLINE void glez_drawstick(real a[], real b[], real r, int nface)
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

#endif
