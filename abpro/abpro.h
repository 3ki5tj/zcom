#include "rv3.h"
#include "rv2.h"

#ifndef ABPRO_H__
#define ABPRO_H__

typedef struct {
  int d; /* dimension */
  int n; /* number of atoms */
  int model; /* 1 or 2 */
  real clj[2][2], sla, slb;
  int *type; /* 0: A, 1: B */
  real *x, *x1, *dx, *dx1;
  real *v;
  real *f;
  real *lmx, *xmin;
  real emin;
} abpro_t;

abpro_t *ab_open(int seqid, int d, int model);
void ab_close(abpro_t *ab);
int ab_checkconn(abpro_t *ab, const real *x);
void ab_shiftcom(abpro_t *ab, real *x);
int ab_writepos(abpro_t *ab, const real *x, const real *v, const char *fname);
int ab_readpos(abpro_t *ab, real *x, real *v, const char *fname);
int ab_initpos(abpro_t *ab);

int ab_shake(abpro_t *ab, const real *x0, real *x1, int itmax, double tol);
int ab_rattle(abpro_t *ab, const real *x0, real *v, int itmax, double tol);
int ab_milcshake(abpro_t *ab, const real *x0, real *x1, real *v, real dt, 
    int itmax, double tol);

real ab_localmin(abpro_t *ab, const real *r, int itmax, double tol);
real ab_energy(abpro_t *ab, const real *r, int soft);
real ab_force(abpro_t *ab, real *f, real *r, int soft);

#endif

