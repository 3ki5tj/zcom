#include "glez.h"

static void display(void)
{
  real a[3] = {-.5f, -.5f, -.5f}, b[3] = {.5f, .5f, .5f}, r = 0.3f;
  float color[4] = {0.8f, 0.6f, 0.0f, 1.0f};
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
  glLightfv(GL_LIGHT0, GL_AMBIENT, color);
  glez_drawstick(a, b, r, 32);
  glutSwapBuffers();
}

static void initgui(void)
{
  GLfloat light_position0[] = {1.0f, 2.0f, 3.0f, 0.0f};  /* Infinite light location. */
  GLfloat bg[] = {0.2f, 0.2f, 0.2f, 1.0f};

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, bg);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
}

int main(int argc, char **argv)
{
  glezInitWindow(&argc, argv, 600, 600, "GLEZ test");
  glutDisplayFunc(display);
  initgui();
  glutMainLoop();
  return 0;
}

