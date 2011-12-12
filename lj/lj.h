#include "def.h"
#include "util.h"
#ifndef LJ_H__
#define LJ_H__

typedef struct {
  int i, j;
  real phi, psi, xi, dx[3], dr2;
} ljpair_t;

typedef struct {
  int d; /* dimension = 3 */
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  real rho;
  real l, vol; /* side length and volume */
  real rc2, rc, rcdef; /* real / preferred rc */

  real * RESTRICT x; /* reduced unit (0, 1) */
  real * RESTRICT v, * RESTRICT f;
  real epot, epots; /* potential energy and shifted potential energy */
  real ekin, tkin, etot;
  real vir, pvir; /* virial and non-ideal gas part of the virial */
  real epot_shift, epot_tail, p_tail;
  double t;

  int usesw; /* switched potential */
  real rs, a4, a5, a6, a7; /* parameters */
  ljpair_t *pr;
  int npr;
  real lap, f2, *gdg, *xdg;

  int usesq; /* square well potential */
  real ra, ra2, rb, rb2; /* -1 for (ra, rb) */
} lj_t;

lj_t *lj_open(int n, int d, real rho, real rcdef);
void lj_close(lj_t *lj);
int lj_writepos(lj_t *lj, const real *x, const real *v, const char *fn);
int lj_readpos(lj_t *lj, real *x, real *v, const char *fn);
real lj_energy(lj_t *lj);
real lj_force(lj_t *lj);
void lj_vv(lj_t *lj, real dt);

#define lj_shiftcom(lj, v)    md_shiftcom(v, lj->n, lj->d)
#define lj_shiftang(lj, x, v) md_shiftang(x, v, lj->n, lj->d)

INLINE void lj_vscale(lj_t *lj, real tp, real ekt) 
 { md_vscale(lj->v, lj->n * lj->d, lj->dof, tp, ekt, &lj->ekin, &lj->tkin); }

INLINE void lj_vrescale(lj_t *lj, real tp, real thermdt) 
 { md_vrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

 /* compute pressure */
INLINE real lj_calcp(lj_t *lj, real tp)
  { return lj->rho * tp + lj->vir / (lj->d * lj->vol) + lj->p_tail; }

void lj_initsw(lj_t *lj, real rs);

/* initialize square well potential */
INLINE void lj_initsq(lj_t *lj, real ra, real rb)
{
  lj->ra = ra;
  lj->ra2 = ra * ra;
  lj->rb = rb;
  lj->rb2 = rb * rb;
  lj->usesq = 1;
  lj_energy(lj);
}

real lj_bconfsw3d(lj_t *lj, real *udb);

#endif

