#ifndef SPECFUNC_C__
#define SPECFUNC_C__

#include "specfunc.h"

/* returns log(Gamma(z)) */
double lngam(double z)
{
  int i;
  double xp, zhg;
  static const double gh = 671./128, sqrt2pi = 2.506628274631000242,
    c[15] = {0.999999999999997092, 57.1562356658629235,-59.5979603554754912,
    14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
    .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
    -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
    .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

  if (z <= 0.) {
    fprintf(stderr, "neg. arg. for gamma(x)\n");
    exit(1);
  }
  for (xp = c[0], i = 1; i < 15; i++)
    xp += c[i]/(z + i);
  zhg = z + gh;
  return (z+.5)*log(zhg) - zhg + log(sqrt2pi*xp/z); /* gamma(z) = gamma(z+1)/z */
}

/* return the p-value, or 1 - cdf(x), for KS distribution */
double ksq(double x)
{
  double y;
  if (x < 0) {
    fprintf(stderr, "neg. arg. for ksq(x)\n");
    exit(1);
  }
  if (x < 1e-15) {
    return 1.;
  } else if (x < 1.18) {
    x = 1.110720734539591525/x;
    y = exp(-x*x);
    return 1. - 2.25675833419102515*x*(y+pow(y,9)+pow(y,25)+pow(y,49));
  } else {
    y = exp(-x*x*2.);
    return 2.*(y - pow(y, 4) + pow(y, 9));
  }
}

#endif

