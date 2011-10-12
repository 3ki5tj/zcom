#ifndef GLEZ_C__
#define GLEZ_C__
#include "glez.h"

/* standard reshape function for GLUT 
 * x: (-1, 1), y: (-1, 1), z: (-10, 10) */
void glez_reshape(int w, int h)
{
  double xs = 1, ys = 1, zs = 10;

  if (w > h) xs = 1.*w/h;
  else ys = 1.*h/w;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-xs, xs, -ys, ys, -zs, zs);
  glMatrixMode(GL_MODELVIEW);

  if (glez_user_reshape) (*glez_user_reshape)(w, h);
}

/* menu function */
void glez_menufunc(int id)
{
  GLfloat mat[4][4];

  if (id <= GLEZ_MENU0 || id >= GLEZ_MENU1) {
    if (glez_user_menufunc) (*glez_user_menufunc)(id);
  } else if (id == GLEZ_SCLM || id == GLEZ_SCLP) {
    GLfloat s = (id == GLEZ_SCLP ? 1.02f : 0.98039216f);
    glScalef(s, s, s);
    glutPostRedisplay();
  } else if (id == GLEZ_ROTXM || id == GLEZ_ROTXP) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == GLEZ_ROTXP) ? 5.f : -5.f, mat[0][0], mat[1][0], mat[2][0]);
    glutPostRedisplay();
  } else if (id == GLEZ_ROTYM || id == GLEZ_ROTYP) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == GLEZ_ROTYP) ? 5.f : -5.f, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  } else if (id == GLEZ_ROTZM || id == GLEZ_ROTZP) {
    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef((id == GLEZ_ROTZP) ? 5.f : -5.f, mat[0][2], mat[1][2], mat[2][2]);
    glutPostRedisplay();
  } else if (id == GLEZ_FULLS) {
    glez_fullscreen();
  } else if (id == GLEZ_QUIT) {
    exit(0);
  }
}

void glez_keyboardfunc(unsigned char c, int x, int y)
{
  int i, j;
  glez_menukey_t *mk;

  (void) x; (void) y;
  if (c == 27) exit(0);
  for (j = 0; j < 2; j++) {
    mk = (j == 0) ? glez_user_menukey : glez_menukey;
    if (mk == NULL) continue;
    for (i = 0; i < mk[i].key; i++)
      if (c == mk[i].key) {
        glez_menufunc(mk[i].id);
        return;
      }
  }
  if (glez_user_keyboardfunc) (*glez_user_keyboardfunc)(c, x, y);
}

void glezMenuKeyFunc(void (*menuf)(int), void (*keyf)(unsigned char, int, int),
    glez_menukey_t *mk)
{
  int i, j;
  char s[256];

  glez_user_menufunc = menuf;
  glez_user_keyboardfunc = keyf;
  glez_user_menukey = mk;
  glutCreateMenu(glez_menufunc);
  for (j = 0; j < 2; j++) {
    for (i = 0; mk[i].key && mk[i].desc; i++) { /* add user menus first */
      sprintf(s, "%s, key: %c\n", mk[i].desc, mk[i].key);
      glutAddMenuEntry(s, mk[i].id);
    }
    mk = glez_menukey; /* switch to system menu */
  }
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutKeyboardFunc(glez_keyboardfunc);
}

void glez_mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    glez_msdown++;
    if (button == 3) {
      glez_menufunc(GLEZ_SCLP);
    } else if (button == 4) {
      glez_menufunc(GLEZ_SCLM);
    }
  } else if (--glez_msdown <= 0) {
    glez_msdown = 0;
  }
  glez_x = x;
  glez_y = y;
  if (glez_user_mouse) (*glez_user_mouse)(button, state, x, y);
}

/* mouse motion function for GLUT */
void glez_motion(int x, int y) 
{
  if (x == glez_x && y == glez_y) return;
  if (glez_msdown) {
    float angx = (float)( (y - glez_y) * 360.f / glutGet(GLUT_WINDOW_HEIGHT) );
    float angy = (float)( (x - glez_x) * 360.f / glutGet(GLUT_WINDOW_WIDTH) );
    float mat[4][4];

    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotated(angx, mat[0][0], mat[1][0], mat[2][0]);
    glRotated(angy, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  }
  glez_x = x;
  glez_y = y;
  if (glez_user_motion) (*glez_user_motion)(x, y);
}

/* toggle full screen state */
void glez_fullscreen(void)
{
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
}

/* draw a stick from a to b, with radius r and nface faces */
void glez_drawstick(real a[], real b[], real r, int nface)
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
