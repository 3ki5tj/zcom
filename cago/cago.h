#include "rv3.h"
#include "util.h"
#include "av.h"
#include "rng.h"
#include "md.h"
#include "pdb.h"
#ifndef CAGO_H__
#define CAGO_H__

/* alpha-carbon based Go-model
 * J. Mol. Biol, Vol. 298 (2000) 937-953 */
typedef struct {
  int n; /* number of residues */
  int dof; /* degree of freedom */
  real kb; /* .5 kb (b - b0)^2, usually 200 */
  real ka; /* .5 ka (a - a0)^2, usually 40 */
  real kd1, kd3; /* kd1 (1 - cos(d - d0)) + kd3 (1 - cos(3*(d-d0))), 1 & 0.5 */
  real nbe, nbc; /* nbc = 4 A */
  unsigned flags; /* input flags */

  rv3_t *xref;
  real epotref; /* energy of the reference structure */
  int *aa;
  real *bref; /* bonds */
  real *aref; /* angle */
  real *dref; /* dihedral */
  real *r2ref; /* pair distance */
  int ncont; /* number of defined contacts */
  int *iscont;

  /* variables for MD simulations */
  rv3_t *x, *v, *f, *x1;
  real ekin, tkin, epot, etot, t;
  real rmsd; /* root-mean-square deviation, result from a rotfit call  */
} cago_t;

#define CAGO_VERBOSE 0x1000
#define CAGO_NCWCA   0x2000  /* use WCA-potential for non-contact pairs,
                              * otherwise, use r^12 repulsion */

#define cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc) \
        cago_openx(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc, \
                  CAGO_VERBOSE + (unsigned) PDB_CONTACT_HEAVY)
/* by default, we use heavy atoms to define contacts,
 * and we are print out the contact information */

cago_t *cago_openx(const char *fnpdb, real kb, real ka, real kd1, real kd3,
    real nbe, real nbc, real rcc, unsigned flags);
void cago_close(cago_t *go);

INLINE void cago_rmcom(cago_t *go, rv3_t *x, rv3_t *v);
INLINE int cago_initmd(cago_t *go, double rndamp, double T0);
INLINE real cago_force(cago_t *go, rv3_t *x, rv3_t *f);
INLINE int cago_vv(cago_t *go, real fscal, real dt);
INLINE real cago_ekin(cago_t *go, rv3_t *v)
  { return go->ekin = md_ekin((real *) v, go->n*3, go->dof, &go->tkin); }
INLINE void cago_vrescale(cago_t *go, real tp, real dt)
  { md_vrescale3d(go->v, go->n, go->dof, tp, dt, &go->ekin, &go->tkin); }
INLINE int cago_mcvrescale(cago_t *go, real tp, real dt)
  { return md_mcvrescale3d(go->v, go->n, go->dof, tp, dt, &go->ekin, &go->tkin); }
INLINE int cago_metro(cago_t *go, real amp, real bet);
INLINE int cago_writepos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn);
INLINE int cago_readpos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn);
INLINE int cago_writepdb(cago_t *go, rv3_t *x, const char *fn);

/* convenient macro for computing RMSD from the reference structure */
#define cago_rmsd(go, x, xf) rv3_rmsd(x, xf, go->xref, NULL, go->n, NULL, NULL)
#define cago_rotfit(go, x, xf) { go->rmsd = cago_rmsd(go, x, xf); }

/* compute the number of contacts from the current configuration */
INLINE int cago_ncontacts(cago_t *go, rv3_t *x, real gam, real *Q, int *mat);

/* copy position or velocities, from `s' to `t' */
#define cago_copyvec(go, t, s) rv3_ncopy(t, s, go->n)

#endif
