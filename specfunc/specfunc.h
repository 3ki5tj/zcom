#ifndef INLINE
#define INLINE static __inline
#endif
#include "util.h"
#ifndef SPECFUNC_H__
#define SPECFUNC_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

INLINE double lngam(double a);
INLINE double lnincgam(double a, double x);
INLINE double lnincgamup(double a, double x);

INLINE double ksq(double x);
INLINE double plegendre(double x, int l, int m);

#endif

