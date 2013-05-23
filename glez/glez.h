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

int glezInitWindow(int *argc, char **argv, int w, int h, const char *name);

static void (*glez_user_reshapefunc)(int w, int h) = NULL;
#define glezReshapeFunc(f)  { glez_user_reshapefunc = f;  glutReshapeFunc(glez_reshapefunc); }

static void (*glez_user_menufunc)(int id) = NULL;
static void (*glez_user_keyboardfunc)(unsigned char c, int x, int y) = NULL;
void glezMenuKeyFunc(void (*menuf)(int), void (*keyf)(unsigned char, int, int),
    glez_menukey_t *mk);

static void (*glez_user_mousefunc)(int button, int state, int w, int h) = NULL;
#define glezMouseFunc(f)    { glez_user_mousefunc = f;  glutMouseFunc(glez_mousefunc); }

static void (*glez_user_motionfunc)(int w, int h) = NULL;
#define glezMotionFunc(f)   { glez_user_motionfunc = f;  glutMotionFunc(glez_motionfunc); }

void glez_fullscreen(void);

void glez_drawstick(real a[], real b[], real r, int nface);

#endif
