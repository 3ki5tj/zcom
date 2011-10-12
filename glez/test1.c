#include "glez.c"

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
  glez_reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
}

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitWindowSize(800, 800);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow("GLEZ test");
  glutDisplayFunc(display);
  glezReshapeFunc(NULL);
  glezMenuKeyFunc(NULL, NULL, NULL);
  glezMouseFunc(NULL);
  glezMotionFunc(NULL);
  initgui();
  glutMainLoop(); 
  return 0;
}

