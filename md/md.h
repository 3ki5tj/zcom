#ifndef ZCINLINE
#define ZCINLINE __inline static
#endif
#include "rv2.h"
#include "rv3.h"
#ifndef MD_H__
#define MD_H__

#define md_shiftcom3d(x, n) md_shiftcom((real *)x, 3, n)
#define md_shiftcom2d(x, n) md_shiftcom((real *)x, 2, n)
#define md_shiftcom(x, d, n) md_shiftcomw(x, NULL, d, n)
#define md_shiftcomw3d(x, w, n) md_shiftcom((real *)x, w, 3, n)
#define md_shiftcomw2d(x, w, n) md_shiftcom((real *)x, w, 2, n)
void md_shiftcomw(real *x, real *w, int d, int n);

void md_shiftang2d(rv2_t *x, rv2_t *v, int n);
void md_shiftang3d(rv3_t *x, rv3_t *v, int n);
ZCINLINE void md_shiftang(real *x, real *v, int d, int n) 
{
  if (d == 2) md_shiftang2d((rv2_t *)x, (rv2_t *)v, n);
  else md_shiftang3d((rv3_t *)x, (rv3_t *)v, n);
}

real md_ekin(real *v, int nd, int dof, real *tkin);

void md_vrescale(real tp, real dt, real *ekin, real *tkin, 
    real *v, int nd, int dof);

#endif

