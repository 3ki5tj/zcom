#ifndef LU_H__
#define LU_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int lusolve(double *a, double *b, int n);
int ludcmp(double *a, int *idx, int n);
int lubksb(double *a, double *b, int *idx, int n);

#endif

