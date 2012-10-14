#ifndef RESTRICT
#define RESTRICT __restrict
#endif
#ifndef INLINE
#define INLINE __inline static
#endif
#include "rv2.h"
#include "rv3.h"
#ifndef MD_H__
#define MD_H__

#define md_shiftcom(x, n, d) md_shiftcomw(x, NULL, n, d)
#define md_shiftcom3d(x, n) md_shiftcomw3d(x, NULL, n)
#define md_shiftcom2d(x, n) md_shiftcomw2d(x, NULL, n)
INLINE void md_shiftcomw(real * RESTRICT x, const real * RESTRICT w, int n, int d);
/* these two are inline instead macros because they offer type checks */
INLINE void md_shiftcomw2d(rv2_t * RESTRICT x, const real * RESTRICT w, int n)
  { md_shiftcomw((real *) x, w, n, 2); }
INLINE void md_shiftcomw3d(rv3_t * RESTRICT x, const real * RESTRICT w, int n)
  { md_shiftcomw((real *) x, w, n, 3); }

INLINE void md_shiftang2d(rv2_t * RESTRICT x, rv2_t * RESTRICT v, int n);
INLINE void md_shiftang3d(rv3_t * RESTRICT x, rv3_t * RESTRICT v, int n);
INLINE void md_shiftang(real * RESTRICT x, real * RESTRICT v, int n, int d) 
{
  if (d == 2) md_shiftang2d((rv2_t *) x, (rv2_t *) v, n);
  else md_shiftang3d((rv3_t *) x, (rv3_t *) v, n);
}

INLINE real md_ekin(const real *v, int nd, int dof, real * RESTRICT tkin);
INLINE real md_ekin2d(rv2_t * RESTRICT v, int n, int dof, real * RESTRICT tkin)
  { return md_ekin((const real *) v, n*2, dof, tkin); }
INLINE real md_ekin3d(rv3_t * RESTRICT v, int n, int dof, real * RESTRICT tkin)
  { return md_ekin((const real *) v, n*3, dof, tkin); }

INLINE void md_vscale(real * RESTRICT v, int nd, int dof, real tp, real ekt, real * RESTRICT ekin, real * RESTRICT tkin);
INLINE void md_vscale2d(rv2_t * RESTRICT v, int n, int dof, real tp, real ekt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vscale((real *) v, n*2, dof, tp, ekt, ekin, tkin); }
INLINE void md_vscale3d(rv3_t * RESTRICT v, int n, int dof, real tp, real ekt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vscale((real *) v, n*3, dof, tp, ekt, ekin, tkin); }

INLINE void md_vrescale(real * RESTRICT v, int nd, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin);
INLINE void md_vrescale2d(rv2_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vrescale((real *) v, n*2, dof, tp, dt, ekin, tkin); }
INLINE void md_vrescale3d(rv3_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vrescale((real *) v, n*3, dof, tp, dt, ekin, tkin); }

INLINE void md_vrescalex(real * RESTRICT v, int nd, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin);
INLINE void md_vrescalex2d(rv2_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vrescalex((real *) v, n*2, dof, tp, dt, ekin, tkin); }
INLINE void md_vrescalex3d(rv3_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { md_vrescalex((real *) v, n*3, dof, tp, dt, ekin, tkin); }

#endif

