#include "rv3.h"
#ifndef CAGO_H__
#define CAGO_H__

/* alpha-carbon based GO-model */
typedef struct {
  int n; /* number of residues */
  int dof; /* degree of freedom */
  real kb; /* .5 kb (b - b0)^2 */
  real ka; /* .5 ka (a - a0)^2 */
  real kd1, kd3; /* kd1 (1 - cos(d - d0)) + kd3 (1 - cos(3*(d-d0))) */
  real nbe, nbc;
  rv3_t *xref;
  int *aa;
  real *bref; /* bonds */
  real *aref; /* angle */
  real *dref; /* dihedral */
  real *r2ref; /* pair distance */
  int *iscont;
  
  /* variables for a typical MD simulation */
  rv3_t *x, *v, *f;
  real ekin, epot, t;
  real rmsd; /* result from a rotfit call  */
} cago_t;

cago_t *cago_open(const char *fnpdb,
  real kb, real ka, real kd1, real kd3, real nbe, real nbc, real rc);
void cago_close(cago_t *go);
void cago_rmcom(cago_t *go, rv3_t *x, rv3_t *v);
int cago_initmd(cago_t *go, double rndamp, double T0);
real cago_force(cago_t *go, rv3_t *x, rv3_t *f);
int cago_vv(cago_t *go, real fscal, real dt);
ZCINLINE real cago_ekin(cago_t *go, rv3_t *v)
  { return md_ekin((real *)v, go->n*3, go->dof, NULL); }
ZCINLINE void cago_vrescale(cago_t *go, real tp, real dt)
  { md_vrescale(tp, dt, &go->ekin, NULL, (real *)go->v, go->n*3, go->dof); }
int cago_writepos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn);
int cago_readpos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn);
int cago_writepdb(cago_t *go, rv3_t *x, const char *fn);

/* convenient macro for computing rmsd from the reference structure */
ZCINLINE real cago_rotfit(cago_t *go, rv3_t *x, rv3_t *xf)
  { return go->rmsd = rotfit3(x, xf, go->xref, NULL, go->n, NULL, NULL); }

#endif
