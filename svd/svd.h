#ifndef SVD_H__
#define SVD_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
int svd(double *a, double *s, double *v, int m, int n);
int svdback(double *u, double *w, double *v, int m, int n, double *x, double *b);
int svdsolve(double *a, double *x, double *b, int n, double tol);

#endif

