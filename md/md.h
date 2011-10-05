#ifndef ZCINLINE
#define ZCINLINE __inline static
#endif
#include "rv2.h"
#include "rv3.h"
#ifndef MD_H__
#define MD_H__

#define md_shiftcom(x, n, d) md_shiftcomw(x, NULL, n, d)
#define md_shiftcom3d(x, n) md_shiftcomw3d(x, NULL, n)
#define md_shiftcom2d(x, n) md_shiftcomw2d(x, NULL, n)
void md_shiftcomw(real *x, const real *w, int n, int d);
/* these two are inline instead macros because they offer type checks */
ZCINLINE void md_shiftcomw2d(rv2_t *x, const real *w, int n) \
  { md_shiftcomw((real *) x, w, n, 2); }
ZCINLINE void md_shiftcomw3d(rv3_t *x, const real *w, int n) \
  { md_shiftcomw((real *) x, w, n, 3); }

void md_shiftang2d(rv2_t *x, rv2_t *v, int n);
void md_shiftang3d(rv3_t *x, rv3_t *v, int n);
ZCINLINE void md_shiftang(real *x, real *v, int n, int d) 
{
  if (d == 2) md_shiftang2d((rv2_t *)x, (rv2_t *)v, n);
  else md_shiftang3d((rv3_t *)x, (rv3_t *)v, n);
}

real md_ekin(const real *v, int nd, int dof, real *tkin);
ZCINLINE real md_ekin2d(rv3_t *v, int n, int dof, real *tkin) \
  { return md_ekin((const real *) v, n*2, dof, tkin); }
ZCINLINE real md_ekin3d(rv3_t *v, int n, int dof, real *tkin) \
  { return md_ekin((const real *) v, n*3, dof, tkin); }

void md_vrescale(real *v, int nd, int dof, real tp, real dt, real *ekin, real *tkin);
ZCINLINE void md_vrescale2d(rv3_t *v, int n, int dof, real tp, real dt, real *ekin, real *tkin) \
    { md_vrescale((real *) v, n*2, dof, tp, dt, ekin, tkin); }
ZCINLINE void md_vrescale3d(rv3_t *v, int n, int dof, real tp, real dt, real *ekin, real *tkin) \
    { md_vrescale((real *) v, n*3, dof, tp, dt, ekin, tkin); }

#endif

