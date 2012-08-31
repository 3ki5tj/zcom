#include "util.h"
#include "av.h"
#include "hist.h"

#ifndef RPT_H__
#define RPT_H__
#ifndef RPT_INF
#define RPT_INF 1e30 /* infinity temperature */
#endif

/* random perturbation temperature, continuous case */
typedef struct  {
  av_t av; /* simple average data */
  hist_t *hs; /* histogram data */
} rpt_t;

INLINE rpt_t *rpt_open(double emin, double emax, double edel)
{
  rpt_t *t;

  xnew(t, 1);
  av_clear(&t->av);
  t->hs = hs_open(1, emin, emax, edel);
  die_if (t->hs == NULL, "failed to initialize rpt\n");
  return t;
}

INLINE void rpt_close(rpt_t *t)
{
  hs_close(t->hs);
  free(t);
}

INLINE void rpt_add(rpt_t *t, double e)
{
  av_add(&(t->av), e);
  hs_add1(t->hs, 0, e, 1.f, HIST_VERBOSE);
}

/* return the temperature */
INLINE double rpt_bet0(const rpt_t *t)
{
  double ave = av_getave(&(t->av)), var = av_getvar(&(t->av));
  var += ave*ave;
  return (var > 0.) ? (2.0*ave/var) : 0.;
}

INLINE double rpt_bet1(const rpt_t *t)
{
  double ave = av_getave(&(t->av)), var = av_getvar(&(t->av));
  return (var > 0.) ? (2.0*ave/var) : 0.;
}

/* obtain a rough estimate of bet, return 1 if it's a special case */
INLINE int rpt_prepbet(const rpt_t *t, double *bet)
{
  int i, l = 0, r = 0, verbose = 0;
  double cnt, sm1, sm2, e, eps = 1e-10, *h = t->hs->arr;
  hist_t *hs = t->hs;

  *bet = 0.;
  /* count the total number and get a rough estimate */
  if (verbose >= 2) printf("  e   #\n");
  for (cnt = sm1 = sm2 = 0, i = 0; i < hs->n; i++) {
    e = hs->xmin + (i + .5) * hs->dx;
    if (h[i] <= 0.) continue;
    if (e < -eps) l = 1; else if (e > eps) r = 1;
    cnt += h[i];
    sm1 += h[i] * e;
    sm2 += h[i] * e * e;
    if (verbose >= 2) printf("%4g %g\n", e, h[i]);
  }
  if (verbose >= 2) printf("\n");
  if (cnt <= 0.) return 1; /* no data */
  /* one-sided, beta = +/- inf. */
  if (!l) { *bet =  RPT_INF; return 1; }
  if (!r) { *bet = -RPT_INF; return 1; }
  if (fabs(sm1) < eps * cnt) return 1; /* even distribution, beta = 0 */
  *bet = 2.0*sm1/sm2; /* should be an underestimate */
  if (verbose)
    printf("Compare: %g, %g;  %g, %g;  %g, %g\n", t->av.s, cnt, t->av.sx, sm1, t->av.sx2, sm2);
  return 0;
}

INLINE double rpt_refinebet(const rpt_t *t, double bet0,
  double (*func)(const hist_t*, double, int, double*), int ord)
{
  int it, verbose = 0;
  double bet, dbet, dbmax, f, df;

  /* increase bet, such that <exp(-bet*de)>  > 1 */
  for (bet = bet0; ; bet *= 1.4)
    if ((f = (*func)(t->hs, bet, ord, &df)) > 0) break;
  if (verbose) printf("expanded bet from %g to %g, f %g, df %g\n", bet0, bet, f, df);

  dbmax = fabs(bet*0.1); /* limit the maximal amount of updating */

  /* iteratively refine the temperature */
  for (it = 0; it <= 1000; it++) {
    f = (*func)(t->hs, bet, ord, &df);
    if (fabs(f) < 1e-14) break;
    dbet = dblconfine(f/df, -dbmax, dbmax);
    if (verbose) printf("f %g, df %g, bet %g, dbet %g\n", f, df, bet, dbet);
    bet += dbet;  /* f should be one, derivative is -df */
  }
  if (verbose >= 2) printf("\n\n");
  return bet;
}

/* evaluate f = < exp(-bet*e) > - 1 and -df/dbet */
INLINE double rpt_getf(const hist_t *hs, double bet, int ord, double *df)
{
  int i;
  double f, cnt, e, xp, *h = hs->arr;

  (void) ord;
  for (f = *df = cnt = 0., i = 0; i < hs->n; i++) {
    if (h[i] <= 0.) continue;
    e = hs->xmin + (i + .5) * hs->dx;
    xp = exp(-bet*e);
    cnt += h[i];
    f   += h[i] * xp;
    *df += h[i] * xp * e;
  }
  f /= cnt;
  *df /= cnt;
  return f - 1;
}

/* obtain the nontrivial solution of < exp(-bet * e) > = 1 */
INLINE double rpt_bet(const rpt_t *t)
{
  double bet;

  if (rpt_prepbet(t, &bet)) return bet;
  return rpt_refinebet(t, bet, rpt_getf, 0);
}

/* evaluate f = (-bet) < min{1, exp(-bet * e)} sgn(e) |e|^ord > and -df/dbet */
INLINE double rpt_getfs(const hist_t *hs, double bet, int ord, double *df)
{
  int i, sgn;
  double f, cnt, e, ep, xp, *h = hs->arr;

  for (f = *df = cnt = 0., i = 0; i < hs->n; i++) {
    if (h[i] <= 0.) continue;
    cnt += h[i];
    e = hs->xmin + (i + .5) * hs->dx;
    sgn = (e > 0.) ? 1 : (e < 0.) ? -1 : 0;
    ep = pow(fabs(e), ord) * sgn;
    xp = -bet*e;
    if (xp >= 0.) {
      xp   = ep;
      f   += h[i] * xp;
    } else {
      xp   = exp(xp) * ep;
      f   += h[i] * xp;
      *df += h[i] * xp * e;
    }
  }
  f *= -bet/cnt;
  *df *= -bet/cnt;
  return f;
}

/* obtain the nontrivial solution of < min{1, exp(-bet * e)} sgn(e) |e|^ord > = 0 */
INLINE double rpt_bets(const rpt_t *t, int ord)
{
  double bet;

  if (rpt_prepbet(t, &bet)) return bet;
  return rpt_refinebet(t, bet, rpt_getfs, ord);
}

/* evaluate f = (-bet) < exp(-bet/2 * e) sgn(e) |e|^ord > */
INLINE double rpt_getfh(const hist_t *hs, double bet, int ord, double *df)
{
  int i, sgn;
  double f, cnt, e, ep, xp, *h = hs->arr;

  for (f = *df = cnt = 0., i = 0; i < hs->n; i++) {
    if (h[i] <= 0.) continue;
    cnt += h[i];
    e = hs->xmin + (i + .5) * hs->dx;
    sgn = (e > 0.) ? 1 : (e < 0.) ? -1 : 0;
    ep = pow(fabs(e), ord) * sgn;
    xp   = exp(-.5*bet*e) * ep;
    f   += h[i] * xp;
    *df += .5* h[i] * xp * e;
  }
  f *= -bet/cnt;
  *df *= -bet/cnt;
  return f;
}

/* obtain the nontrivial solution of < exp(-bet/2 * e) sgn(e) |e|^ord > = 0 */
INLINE double rpt_beth(const rpt_t *t, int ord)
{
  double bet;

  if (rpt_prepbet(t, &bet)) return bet;
  return rpt_refinebet(t, bet, rpt_getfh, ord);
}

/* write distribution to file */
INLINE int rpt_wdist(rpt_t *t, const char *fn)
{
  return hs_save(t->hs, fn, HIST_ADDAHALF);
}

#define RPTI_HASINF 0x1000

typedef struct {
  int emin, emax, edel, m;
  int *h; /* histogram */
  unsigned flags;
} rpti_t; /* perturbation temperature, integer */

INLINE rpti_t *rpti_open(int emin, int emax, int edel, unsigned flags)
{
  rpti_t *t;

  xnew(t, 1);
  die_if (emin >= emax, "emin %d >= emax %d\n", emin, emax);
  t->emin = emin;
  t->emax = emax;
  t->edel = edel;
  t->m = (emax - emin)/edel + 1;
  t->flags = flags;
  xnew(t->h, t->m + 1); /* add one for infinity */
  return t;
}

INLINE void rpti_close(rpti_t *t)
{
  free(t->h);
  free(t);
}

/* add an energy increment */
INLINE void rpti_add(rpti_t *t, int e)
{
  int ie;

  die_if (e < t->emin, "e %d smaller than emin %d\n", e, t->emin);
  ie = (e - t->emin) / t->edel;
  if (e > t->emax) {
    if (t->flags & RPTI_HASINF)
      t->h[t->m]++;
    else
      die_if (e > t->emax, "e %d larger than emax %d\n", e, t->emax);
  } else {
    t->h[ie]++;
  }
}


/* count data points */
INLINE double rpti_cnt(rpti_t *t)
{
  int i;
  double cnt = 0;
  for (i = 0; i < t->m + 1;  i++)
    cnt += t->h[i];
  return cnt;
}


/* estimate the temperature from 2 <e> / <e^2>, and 2 <e> / <De^2> */
INLINE double rpti_bet1(const rpti_t *t, double *bet0)
{
  int i;
  double e, cnt = 0, sm = 0., sm2 = 0.;

  for (i = 0; i < t->m;  i++) { /* do not count infinity at i == t->m */
    e = t->emin + i * t->edel; /* the energy change */
    cnt += t->h[i];
    sm  += t->h[i] * e;
    sm2 += t->h[i] * e * e;
  }
  if (bet0) *bet0 = 0.;
  if (cnt <= 0.) return 0.;
  sm /= cnt;
  sm2 = sm2/cnt;
  if (t->h[t->m] > 0 && sm < 0.) sm = 0.; /* has infinity */
  if (bet0) *bet0 = 2.0*sm/sm2;
  sm2 -= sm*sm;
  if (sm2 <= 0.) return (sm > 0.) ? RPT_INF : -RPT_INF;
  return 2.0*sm/sm2;
}

/* a rough estimate for bet, return 1 if it's a special case  */
INLINE rpti_prepbet(const rpti_t *t, double *bet)
{
  int i, e;
  int hl, hr;

  /* I. get a rough estimate */
  rpti_bet1(t, bet);
  if (fabs(*bet) < 1e-14 || fabs(*bet) > .99*RPT_INF) return 1;

  /* II. check if the distribution is single-sided, which means infinite beta */
  for (hl = hr = 0, i = 0; i <= t->m; i++) {
    e = t->emin + i * t->edel;
    if (e < 0) hl += t->h[i];
    else if (e > 0) hr += t->h[i];
  }
  if (hl <= 0.) { *bet = RPT_INF; return 1; }
  if (hr <= 0.) { *bet =-RPT_INF; return 1; }
  return 0;
}

/* adjust beta until (*func) becomes 0 */
INLINE double rpti_refinebet(const rpti_t *t, double bet0,
  double (*func)(const rpti_t*, double, int, double*), int ord)
{
  double bet, dbmax, dbet, f, df;
  int it, verbose = 0;

  /* I. increase |bet|, such that <func> > 0 */
  for (bet = bet0; ; bet *= 1.4)
    if ((f = (*func)(t, bet, ord, &df)) > 0) break;
  if (verbose) printf("expanded bet from %g to %g, f %g, df %g\n", bet0, bet, f, df);

  dbmax = fabs(bet*0.1); /* limit the maximal amount of updating */

  /* II. iteratively refine the temperature */
  for (it = 0; it <= 1000; it++) {
    f = (*func)(t, bet, ord, &df);
    /* f should be one, derivative is -df */
    if (fabs(f) < 1e-14) break;
    dbet = dblconfine(f/df, -dbmax, dbmax);
    if (verbose) printf("f %g, df %g, bet %g, dbet %g\n", f, df, bet, dbet);
    bet += dbet;
  }
  if (verbose) printf("\n\n");
  return bet;
}

/* evaluate f = < exp(-bet*e) - 1> and -df/dbet */
INLINE double rpti_getf(const rpti_t *t, double bet, int ord, double *df)
{
  int i, e;
  double f, cnt, xp;

  (void) ord; /* used */
  for (f = *df = cnt = 0., i = 0; i < t->m; i++) {
    e = t->emin + i * t->edel;
    if (t->h[i] <= 0) continue;
    xp = exp(-bet*e);
    cnt += t->h[i];
    f += t->h[i] * xp;
    *df += t->h[i] * xp * e;
  }
  /* we should not include the infinity here */
  f /= cnt;
  *df /= cnt;
  return f - 1;
}

/* estimated the temperature using the identity approach */
INLINE double rpti_bet(rpti_t *t)
{
  double bet;

  if (rpti_prepbet(t, &bet)) return bet;
  return rpti_refinebet(t, bet, rpti_getf, 0);
}

/* evaluate f = (-bet) < min{1, exp(-bet*e)} sgn(e) * |e|^ord > and -d f/d bet */
INLINE double rpti_getfs(const rpti_t *t, double bet, int ord, double *df)
{
  int i, e, sgn, k, ep;
  double f, cnt, xp;

  die_if (bet < 0. && t->h[t->m] > 0, "bet %g must be positive with inf. (%d)\n", bet, t->h[t->m]);
  for (f = *df = cnt = 0., i = 0; i < t->m; i++) {
    if (t->h[i] <= 0.) continue;
    cnt += t->h[i];
    e = t->emin + i * t->edel;
    sgn = (e > 0) ? 1 : (e < 0) ? -1 : 0;
    /*  double ep = pow(fabs(e), ord) * sgn; */
    for (ep = 1, k = 0; k < ord; k++) ep *= e;
    ep = abs(ep) * sgn;
    xp = -bet*e;
    if (xp >= 0.) {
      xp     = ep;
      f     += t->h[i] * ep;
    } else {
      xp     = exp(xp) * ep;
      f     += t->h[i] * xp;
      *df   += t->h[i] * xp * e;
    }
  }
  f *= -bet/cnt;
  *df *= -bet/cnt;
  return f;
}

/* estimated the temperature using area integrals */
INLINE double rpti_bets(rpti_t *t, int ord)
{
  double bet;

  if (rpti_prepbet(t, &bet)) return bet;
  return rpti_refinebet(t, bet, rpti_getfs, ord);
}


/* evaluate f = (-bet) < exp(-bet/2*e) sgn(e) * |e|^ord > and -d f/d bet */
INLINE double rpti_getfh(const rpti_t *t, double bet, int ord, double *df)
{
  int i, e, sgn, k, ep;
  double f, cnt, xp;

  die_if (bet < 0. && t->h[t->m] > 0, "bet %g must be positive with inf. (%d)\n", bet, t->h[t->m]);
  for (f = *df = cnt = 0., i = 0; i < t->m; i++) {
    if (t->h[i] <= 0.) continue;
    cnt += t->h[i];
    e = t->emin + i * t->edel;
    sgn = (e > 0) ? 1 : (e < 0) ? -1 : 0;
    /*  double ep = pow(fabs(e), ord) * sgn; */
    for (ep = 1, k = 0; k < ord; k++) ep *= e;
    ep = abs(ep) * sgn;
    xp     = exp(-.5*bet*e) * ep;
    f     += t->h[i] * xp;
    *df   += .5 * t->h[i] * xp * e;
  }
  f *= -bet/cnt;
  *df *= -bet/cnt;
  return f;
}

/* estimated the temperature using area integrals */
INLINE double rpti_beth(rpti_t *t, int ord)
{
  double bet;

  if (rpti_prepbet(t, &bet)) return bet;
  return rpti_refinebet(t, bet, rpti_getfh, ord);
}


/* write distribution to file */
INLINE int rpti_wdist(rpti_t *t, const char *fn)
{
  int i, i0 = -1, i1 = -1, cnt = 0;
  FILE *fp;

  for (i = 0; i < t->m; i++) {
    if (t->h[i] > 0) {
      if (i0 < 0) i0 = i;
      if (i > i1) i1 = i;
      cnt += t->h[i];
    }
  }
  /* cnt += t->h[t->m]; */

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# 0 0 | 1 %d %d %d |\n", t->m, t->emin, t->edel);
  for (i = i0; i <= i1; i++)
    fprintf(fp, "%d %d %.8f 0\n",
      t->emin + i * t->edel, t->h[i], 1.0*t->h[i]/t->edel/cnt);
  if (t->h[i = t->m] > 0)
    fprintf(fp, "inf %d %.8f 0\n", t->h[i], 1.0*t->h[i]/t->edel/cnt);
  fclose(fp);
  return 0;
}

#endif /* RPT_H__ */

