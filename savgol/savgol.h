#ifndef SAVGOL_H__
#define SAVGOL_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *savgol(int w, int ord, int der, int h, int verbose);
double *savgol2d(int iw, int jw, int ord, int h, int verbose);

#endif

