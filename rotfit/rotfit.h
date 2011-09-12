#include "rv3.h"
#ifndef ROTFIT_H__
#define ROTFIT_H__

real rotfit3(rv3_t *x, rv3_t *xf, rv3_t *y, const real *w, int n, 
    real (*r)[3], real *t);

#endif
