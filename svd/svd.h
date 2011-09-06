#include "def.h"
#ifndef SVD_H__
#define SVD_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
int svd(real *a, real *s, real *v, int m, int n);
int svdback(real *u, real *w, real *v, int m, int n, real *x, real *b);
int svdsolve(real *a, real *x, real *b, int n, real tol);

#endif

