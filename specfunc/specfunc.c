#ifndef SPECFUNC_C__
#define SPECFUNC_C__

#include "specfunc.h"

/* returns log(Gamma(a)),
 * where Gamma(a) = \int_0^\infty e^(-t) t^(a-1) dt */
INLINE double lngam(double a)
{
  int i;
  double xp, ahg;
  static const double gh = 671./128, sqrt2pi = 2.506628274631000242,
    c[15] = {0.999999999999997092, 57.1562356658629235,-59.5979603554754912,
    14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
    .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
    -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
    .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

  die_if (a <= 0., "neg. arg. for lngam(%g)\n", a);
  for (xp = c[0], i = 1; i < 15; i++)
    xp += c[i]/(a + i);
  ahg = a + gh;
  return (a+.5)*log(ahg) - ahg + log(sqrt2pi*xp/a); /* gamma(a) = gamma(a+1)/a */
}

/* returns incomplete gamma function log(gamma(a, x)),
 * where gamma(a, x) = \int_0^x e^(-t) t^(a-1) dt
 * = e^(-x) x^a \sum_{i = 0} Gamma(a)/Gamma(a + 1 + i) x^i
 * for small x, cf. gser() in Numerical-Recipes */
INLINE double lnincgam0(double a, double x)
{
  int i;
  double del, sum;

  die_if (x < 0., "neg. arg. for lnincgam0(%g)\n", x);
  if (x <= 0) return -1e30; /* log(0+) */
  sum = del = 1.0/a;
  for (i = 1; i <= 1000; i++) {
    del *= x/(a + i);
    sum += del;
    if (fabs(del) < fabs(sum)*5e-16) break;
  }
  return -x + a * log(x) + log(sum);
}

/* returns incomplete gamma function log(Gamma(a, x)),
 * where Gamma(a, x) = \int_x^\infty e^(-t) t^(a-1) dt
 * = e^(-x) x^a [1/(x+1-a-) [1*(1-a)/(x+3-a-) [2*(2-a)/(x+5-a) ...]
 * for large x, cf. gcf() in Numerical-Recipes */
INLINE double lnincgam1(double a, double x)
{
  int i;
  double an, b, c, d, h, del;
  const double fpmin = 1e-300;

  b = x + 1 - a;
  d = 1/b;
  c = 1/fpmin;
  h = d;
  /* modified Lentz's method for the continued fraction */
  for (i = 1; i <= 1000; i++) {
    an = -1.*i*(i - a); /* numerator */
    b += 2;
    d = d*an + b;
    if (fabs(d) < fpmin) d = fpmin;
    c = b + an/c;
    if (fabs(c) < fpmin) c = fpmin;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del - 1) < 5e-16) break;
  }
  return -x + a * log(x) + log(h);
}

/* returns incomplete gamma function log(gamma(a, x)),
 * where gamma(a, x) = \int_0^x e^(-t) t^(a-1) dt */
INLINE double lnincgam(double a, double x)
{
  return (x < a + 1) ? lnincgam0(a, x) : lndif(lngam(a), lnincgam1(a, x));
}

/* returns incomplete gamma function log(Gamma(a, x)),
 * where Gamma(a, x) = \int_x^\infty e^(-t) t^(a-1) dt */
INLINE double lnincgamup(double a, double x)
{
  return (x < a + 1) ? lndif(lngam(a), lnincgam0(a, x)) : lnincgam1(a, x);
}

/* return the p-value, or 1 - cdf(x), for KS distribution */
INLINE double ksq(double x)
{
  double y;
  die_if (x < 0, "neg. arg. for ksq(x = %g)\n", x);
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

/* normalized associated legendre polynomial
 * nplm(x) = sqrt( (2l+1)/(4 pi) (l-m)!/(l+m)!) P_l^m(x)
 * real solution of m <= l, l, m >= 0
 * (1 - x^2) y'' - 2 x y' + [l(l+1) - m^2/(1-x^2)] y = 0 */
INLINE double plegendre(double x, int l, int m)
{
  int i;
  double y, yp, ypp, f, fp, s = 1 - x*x;
 
  if (m < 0 || m > l || s < 0) return 0; 
  for (yp = 1, i = 1; i <= m; i++) yp *= (1 + .5/i)*s;
  yp = sqrt(yp/(4*M_PI)) * (m % 2 ? -1: 1); /* P(m, m) */
  /* (l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1)*P_{l-2}^m */
  for (fp = 1, ypp = 0, i = m + 1; i <= l; i++, fp = f, ypp = yp, yp = y) {
    f = sqrt((4.*i*i-1)/((i-m)*(i+m)));
    y = f*(x*yp - ypp/fp);
  }
  return yp;
}

#endif

