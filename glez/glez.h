#include "rv3.h"
#ifndef GLEZ_H__
#define GLEZ_H__
#include <GL/glut.h>

int glez_x, glez_y; /* current position */
int glez_msdown; /* mouse state */
float glez_zoomscale = 1.0f;

typedef struct { 
  int id;
  unsigned char key;
  const char *desc;
} glez_menukey_t;

glez_menukey_t *glez_user_menukey;

enum { GLEZ_MENU0 = 1000, 
  GLEZ_ROTXM, GLEZ_ROTXP, GLEZ_ROTYM, GLEZ_ROTYP, GLEZ_ROTZM, GLEZ_ROTZP, 
  GLEZ_SCLM, GLEZ_SCLP, GLEZ_FULLS, GLEZ_QUIT, GLEZ_MENU1 };

/* common menu-key pairs */
glez_menukey_t glez_menukey[] = {
  {GLEZ_ROTXM,  'x', "Rotate around -x"},
  {GLEZ_ROTXP,  'X', "Rotate around +x"},
  {GLEZ_ROTYM,  'y', "Rotate around -y"},
  {GLEZ_ROTYP,  'Y', "Rotate around +y"},
  {GLEZ_ROTZM,  'z', "Rotate around -z"},
  {GLEZ_ROTZP,  'Z', "Rotate around +z"},
  {GLEZ_SCLM,   's', "Zoom out"},
  {GLEZ_SCLP,   'S', "Zoom in"},
  {GLEZ_FULLS,  'f', "Toggle fullscreen"},
  {GLEZ_QUIT,   'q', "Quit"},
  {-1, '\0', NULL}};

enum {GLEZ_MSNONE, GLEZ_MSZOOM, GLEZ_MSROTATE};

void (*glez_user_reshape)(int w, int h);
void glez_reshape(int w, int h);
#define glezReshapeFunc(f)  { glez_user_reshape = f;  glutReshapeFunc(glez_reshape); }

void (*glez_user_menufunc)(int id);
void glez_menufunc(int id);
void (*glez_user_keyboardfunc)(unsigned char c, int x, int y);
void glez_keyboardfunc(unsigned char c, int x, int y);
void glezMenuKeyFunc(void (*menuf)(int), void (*keyf)(unsigned char, int, int),
    glez_menukey_t *mk);

void (*glez_user_mouse)(int button, int state, int w, int h);
void glez_mouse(int button, int state, int x, int y);
#define glezMouseFunc(f)    { glez_user_mouse   = f;  glutMouseFunc(glez_mouse); }

void (*glez_user_motion)(int w, int h);
void glez_motion(int x, int y); 
#define glezMotionFunc(f)   { glez_user_motion  = f;  glutMotionFunc(glez_motion); }

void glez_fullscreen(void);

void glez_drawstick(real a[], real b[], real r, int nface);

#endif
