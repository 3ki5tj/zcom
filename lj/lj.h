#include "def.h"
#ifndef LJ_H__
#define LJ_H__

typedef struct {
  int d; /* dimension = 3 */
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  real rho;
  real l, vol; /* side length and volume */
  real rc, rcdef; /* real / preferred rc */

  real * RESTRICT x; /* reduced unit (0, 1) */
  real * RESTRICT v, * RESTRICT f;
  real epot, epots; /* potential energy and shifted potential energy */
  real ekin, tkin, etot;
  real vir, pvir; /* virial and non-ideal gas part of the virial */
  real epot_shift, epot_tail, p_tail;
  double t;
} lj_t;

lj_t *lj_open(int n, int d, real rho, real rcdef);
void lj_close(lj_t *lj);
void lj_energy(lj_t *lj);
void lj_force(lj_t *lj);
void lj_vv(lj_t *lj, real dt);

INLINE void lj_vrescale(lj_t *lj, real tp, real thermdt) 
 { md_vrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

 /* compute pressure */
INLINE real lj_calcp(lj_t *lj, real tp)
  { return lj->rho * tp + lj->vir / (lj->d * lj->vol) + lj->p_tail; }

#endif

