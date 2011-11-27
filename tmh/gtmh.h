/* a helper header to display gtmh progress */
#ifndef GTMH_H__
#define GTMH_H__
#if defined(Macintosh) || defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int gtmh_winid; /* window id */
int gtmh_w, gtmh_h; /* window size */
const char *gtmh_name;
int gtmh_pause = 0;
float gtmh_aspect = 1.0f; /* y : x */
float gtmh_xmargin = 0.1f, gtmh_ymargin = 0.05f;
float gtmh_eheight = 0.3f; /* energy histogram height */
float gtmh_dotsz = 0.001f; /* size of a dot */
float gtmh_dhdescale = 2.0f; /* vertical scale for dhde curve */
tmh_t *gtmh_tmh = NULL;
double gtmh_epot = 0.0;
double *gtmh_ehis;

/* colors */
GLfloat gtmh_cdhde[] = {1.0f, 0.0f, 0.0f, 0.f};  /* dhde curve */
GLfloat gtmh_cdhdegrid[] = {0.3f, 0.0f, 0.0f, 0.f}; /* dhde grids */
GLfloat gtmh_cdark[] = {0.1f, 0.3f, 0.5f, 0.f}; /* background energy interval bands */
GLfloat gtmh_cdarker[] = {0.0f, 0.05f, 0.15f, 0.f};
GLfloat gtmh_cgray[] = {0.6f, 0.6f, 0.6f, 0.f};
GLfloat gtmh_ccurr[] = {0.1f, 1.0f, 0.2f, 0.f};
GLfloat gtmh_cehis1[] = {0.5f, 0.5f, 0.0f, 0.f}; /* energy histogram bright color */
GLfloat gtmh_cehis2[] = {0.3f, 0.1f, 0.0f, 0.f}; /* energy histogram dark color */

/* reshape function */
static void gtmh_reshape(int w, int h)
{
  float asp = 1.0f * h / w; /* screen aspect ratio */

  gtmh_w = w;
  gtmh_h = h;
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gtmh_aspect = asp;
  glOrtho(-gtmh_xmargin, 1.f + gtmh_xmargin, 
    -gtmh_ymargin * asp, (1 + gtmh_eheight + gtmh_ymargin)*asp, -5.0f, 5.0f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/* set the current drawing color */
INLINE void gtmh_setcolor(GLfloat c[]) { glColor3fv(c); }

/* convert dhde value to vertical screen scale */
INLINE float gtmh_dhde2y(double x)
{ return (((float) x - 1.0f) * gtmh_dhdescale + 0.5f) * gtmh_aspect; }

/* draw a square dot (for a footprint) */
INLINE void gtmh_drawdot(float x, float y)
{
  y *= gtmh_aspect;
  glBegin(GL_POLYGON);
  glVertex2f(x - gtmh_dotsz, y - gtmh_dotsz);
  glVertex2f(x - gtmh_dotsz, y + gtmh_dotsz);
  glVertex2f(x + gtmh_dotsz, y + gtmh_dotsz);
  glVertex2f(x + gtmh_dotsz, y - gtmh_dotsz);
  glEnd();
}

/* draw the current position and some past point */
static void gtmh_drawfootprint(float x, float y, float w, float h)
{
#define GTMH_XYCNT 1024
  static float xy[GTMH_XYCNT][2]; /* history of trace */
  static int xyfull = 0, xyid = 0;
  float c[3], f = 0.3f;
  int j, i, ic, icnt;

  /* add the current point to the queue */
  xy[xyid][0] = x;
  xy[xyid][1] = y;
  xyid = (xyid + 1) % GTMH_XYCNT;
  if (xyid == 0) xyfull = 1;

  /* draw past points */
  icnt = xyfull ? GTMH_XYCNT : xyid;
  for (ic = icnt - 1; ic >= 0; ic--) { /* ic is the distance from the current point */
    i = (xyid - ic + GTMH_XYCNT) % GTMH_XYCNT; /* actual index */
    for (j = 0; j < 3; j++) /* past points dimmer */
      c[j] = gtmh_ccurr[j]*(1.f - 1.f*ic/icnt);
    gtmh_setcolor(c);
    gtmh_drawdot(xy[i][0], xy[i][1]);
  }

  /* draw the current point */
  y *= gtmh_aspect;
  gtmh_setcolor(gtmh_ccurr);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  glVertex2f(x - w, y); glVertex2f(x + w, y);
  glVertex2f(x, y - h); glVertex2f(x, y + h);
  glVertex2f(x - w*f, y - h); glVertex2f(x + w*f, y - h);
  glVertex2f(x - w*f, y + h); glVertex2f(x + w*f, y + h);
  glVertex2f(x - w, y - h*f); glVertex2f(x - w, y + h*f);
  glVertex2f(x + w, y - h*f); glVertex2f(x + w, y + h*f);
  glEnd();
}

/* draw energy histogram */
static void gtmh_drawehis(void)
{
  int i;
  double hmax = 0., x, y, y0 = gtmh_aspect, ie;
  tmh_t *m = gtmh_tmh;

  if (gtmh_ehis == NULL) return;
  for (i = 0; i < m->en; i++)
    if (gtmh_ehis[i] > hmax) hmax = gtmh_ehis[i];
  gtmh_setcolor(gtmh_cehis2);
  for (i = 0; i < m->en; i++) {
    glBegin(GL_POLYGON);
    ie = m->emin + i * m->de;
    x = (ie - m->erg0)/(m->erg1 - m->erg0);
    y = gtmh_eheight * gtmh_ehis[i] / hmax;
    glVertex2d(x, y0);
    gtmh_setcolor(gtmh_cehis1);
    glVertex2d(x, y0 + gtmh_aspect * y);
    x = (ie + m->de - m->erg0)/(m->erg1 - m->erg0);
    glVertex2d(x, y0 + gtmh_aspect * y);
    gtmh_setcolor(gtmh_cehis2);
    glVertex2d(x, y0);
    glEnd();
  }
}

/* display */
static void gtmh_display(void)
{
  int i;
  float x, y;
  tmh_t *m = gtmh_tmh;

  if (gtmh_pause) return;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  /* energy background strips */
  for (i = -1; i <= m->ergn; i++) {
    gtmh_setcolor(i % 2 == 0 ? gtmh_cdark : gtmh_cdarker);
    glBegin(GL_POLYGON);
    x = (i >= 0) ? (1.0f*i/m->ergn) : -1.0f;
    glVertex2f(x, 0);
    glVertex2f(x, gtmh_aspect);
    x = (i < m->ergn) ? (1.0f*(i+1)/m->ergn) : 2.0f;
    glVertex2f(x, gtmh_aspect);
    glVertex2f(x, 0);
    glEnd();
  }

  /* draw energy histogram */
  gtmh_drawehis();

  /* axes */
  gtmh_setcolor(gtmh_cgray);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  glVertex2f(0, 0); glVertex2f(1, 0); /* x-axis */
  glVertex2f(0, gtmh_aspect); glVertex2f(1, gtmh_aspect);
  glVertex2f(0, 0); glVertex2f(0, gtmh_aspect); /* y-axis */
  glVertex2f(1, 0); glVertex2f(1, gtmh_aspect);
  glEnd();

  /* reference E-T */
  gtmh_setcolor(gtmh_ccurr);
  glLineWidth(1.0f);
  glBegin(GL_LINES);
  glVertex2f(0.f, 0.f); glVertex2f(1.f, gtmh_aspect);
  glEnd();

  /* dhde grid */
  gtmh_setcolor(gtmh_cdhdegrid);
  glLineStipple(1, 0x5555);
  glEnable(GL_LINE_STIPPLE);
  glBegin(GL_LINES);
  for (x = 0; x < 10.f; x += 0.1f) {
    y = gtmh_dhde2y(x);
    if (y < 0 || y > gtmh_aspect) continue;
    glVertex2f(0, y);
    glVertex2f(1.f, y);
  }
  glEnd();
  glDisable(GL_LINE_STIPPLE);

  /* dhde reference 1.0 */
  gtmh_setcolor(gtmh_cdhde);
  glLineWidth(1.0f);
  glBegin(GL_LINES);
  glVertex2f(0, 0.5f * gtmh_aspect); glVertex2f(1.f, .5f * gtmh_aspect);
  glEnd();

  /* dhde */
  glLineWidth(2.0);
  glBegin(GL_LINE_STRIP);
  for (i = 0; i <= m->ergn; i++) {
    x = 1.f*i/m->ergn;
    y = gtmh_dhde2y(m->dhde[i]);
    glVertex2f(x, y);
  }
  glEnd();
  glLineWidth(1.0);

  /* draw trace */
  x = (float)( (gtmh_epot - m->erg0)/(m->erg1 - m->erg0) );
  y = (float)( (m->tp - m->tp0)/(m->tp1 - m->tp0) );
  gtmh_drawfootprint(x, y, 0.02f * gtmh_aspect, 0.03f * gtmh_aspect);

  glutSwapBuffers();
}

/* create a window */
INLINE void gtmh_init(tmh_t *m, int w, int h, const char *caption)
{
  gtmh_tmh = m;
  xnew(gtmh_ehis, m->en + 1);
  gtmh_w = w;
  gtmh_h = h;
  gtmh_name = caption;
  glutInitWindowSize(w, h);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  gtmh_winid = glutCreateWindow(caption);
  glutReshapeFunc(gtmh_reshape);
  glutDisplayFunc(gtmh_display);
}

/* register the recent epot */
INLINE void gtmh_eadd(double epot)
{
  tmh_t *m = gtmh_tmh;
  if (epot > m->emin && epot < m->emax) {
    int ie = (int)((epot - m->emin)/m->de);
    gtmh_ehis[ie]++;
  }
  gtmh_epot = epot;
}

#endif
