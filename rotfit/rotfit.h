#include "util.h"
#include "rv3.h"
#ifndef ROTFIT_H__
#define ROTFIT_H__

int rotfit(rv3_t *x0, rv3_t *y0, const real *w, int n, real r[3][3], real t[3]);
real rmsdev(rv3_t *x0, rv3_t *x, rv3_t *y0, const real *w, int n, 
    real r[3][3], const real t[3]);

#endif
