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

typedef struct {
  int emin, emax, edel, m;
  int *h; /* histogram */
} rpti_t; /* perturbation temperature, integer */

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

/* evaluate f = < exp(-bet*e) > and -df/dbet */
INLINE double rpt_getf(const hist_t *hs, double bet, double *df)
{
  int i;
  double f, cnt, e, xp, *h = hs->arr;

  for (f = *df = cnt = 0., i = 0; i < hs->n; i++) {
    e = hs->xmin + (i + .5) * hs->dx;
    xp = exp(-bet*e);
    cnt += h[i];
    f   += h[i] * xp;
    *df += h[i] * xp * e;
  }
  f /= cnt;
  *df /= cnt;
  return f;
}

/* obtain the nontrivial solution of < exp(-bet * e) > = 1 */
INLINE double rpt_bet(const rpt_t *t)
{
  int it, i, l = 0, r = 0;
  double bet, dbet, dbmax, cnt, f, df, e, eps = 1e-10, *h = t->hs->arr;
  hist_t *hs = t->hs;
  int verbose = 1;

  /* count the total number and get a rough estimate */
  if (verbose >= 2) printf("  e   #\n");
  for (cnt = f = df = 0, i = 0; i < hs->n; i++) {
    e = hs->xmin + (i + .5) * hs->dx;
    if (h[i] <= 0.) continue;
    if (e < -eps) l = 1; else if (e > eps) r = 1;
    cnt += h[i];
    f   += h[i] * e;
    df  += h[i] * e * e;
    if (verbose >= 2) printf("%4g %g\n", e, h[i]);
  }
  if (verbose >= 2) printf("\n");
  if (cnt <= 0.) return 0.; /* no data */
  /* one-sided, beta = +/- inf. */
  if (!l) return RPT_INF;
  if (!r) return -RPT_INF;
  if (fabs(f) < eps * cnt) return 0.0; /* even distribution, beta = 0 */
  bet = dbet = 2.0*f/df; /* should be an underestimate */

  /* increase bet, such that <exp(-bet*de)>  > 1 */
  for (; ; bet *= 1.4)
    if ((f = rpt_getf(t->hs, bet, &df)) > 1) break;

  if (verbose) {
    printf("Compare: %g, %g;  %g, %g;  %g, %g\n", t->av.s, cnt, t->av.sx, f, t->av.sx2, df);
    printf("first estimate being %g, cf. %g, expanded to %g\n", dbet, rpt_bet0(t), bet);
  }

  dbmax = fabs(bet*0.1); /* limit the maximal amount of updating */

  /* iteratively refine the temperature */
  for (it = 0; it <= 1000; it++) {
    f = rpt_getf(t->hs, bet, &df);
    if (fabs(f - 1) < 1e-14 || fabs(df) < 1e-12) break;
    dbet = dblconfine((f - 1)/df, -dbmax, dbmax);
    if (verbose) printf("f %g, df %g, bet %g, dbet %g\n", f, df, bet, dbet);
    bet += dbet;  /* f should be one, derivative is -df */
  }
  if (verbose >= 2) printf("\n\n");
  return bet;
}

/* write distribution to file */
INLINE int rpt_wdist(rpt_t *t, const char *fn)
{
  return hs_save(t->hs, fn, HIST_ADDAHALF);
}


INLINE rpti_t *rpti_open(int emin, int emax, int edel)
{
  rpti_t *t;

  xnew(t, 1);
  die_if(emin >= emax, "emin %d >= emax %d\n", emin, emax);
  t->emin = emin;
  t->emax = emax;
  t->edel = edel;
  t->m = (emax - emin)/edel + 1;
  xnew(t->h, t->m);
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
  die_if(e < t->emin || e > t->emax, "e %d out of range (%d, %d)\n", e, t->emin, t->emax);
  t->h[(e - t->emin)/t->edel]++;
}


/* count data points */
INLINE double rpti_cnt(rpti_t *t)
{
  int i;
  double cnt = 0;
  for (i = 0; i < t->m;  i++)
    cnt += t->h[i];
  return cnt;
}


/* estimate the temperature from 2 <e> / <e^2>, and 2 <e> / <De^2> */
INLINE double rpti_bet1(rpti_t *t, double *bet0)
{
  int i;
  double e, cnt = 0, sm = 0., sm2 = 0.;
  for (i = 0; i < t->m;  i++) {
    e = t->emin + i * t->edel; /* the energy change */
    cnt += t->h[i];
    sm  += t->h[i] * e;
    sm2 += t->h[i] * e * e;
  }
  if (bet0) *bet0 = 0.;
  if (cnt <= 0.) return 0.;
  sm /= cnt;
  sm2 = sm2/cnt;
  if (bet0) *bet0 = 2.0*sm/sm2;
  sm2 -= sm*sm;
  if (sm2 <= 0.) return (sm > 0.) ? RPT_INF : -RPT_INF;
  return 2.0*sm/sm2;
}

/* evaluate f = < exp(-bet*e) > and -df/dbet */
INLINE double rpti_getf(const rpti_t *t, double bet, double *df)
{
  int i;
  double f, cnt, e, xp;

  for (f = *df = cnt = 0., i = 0; i < t->m; i++) {
    e = t->emin + i * t->edel;
    xp = exp(-bet*e);
    cnt += t->h[i];
    f += t->h[i] * xp;
    *df += t->h[i] * xp * e;
  }
  f /= cnt;
  *df /= cnt;
  return f;
}

/* estimated the temperature using the identity approach */
INLINE double rpti_bet(rpti_t *t)
{
  double bet, dbet, dbmax, f, df;
  int i, e, it, hl, hr, verbose = 0;

  /* I. get a rough estimate */
  bet = rpti_bet1(t, NULL);
  if (fabs(bet) < 1e-14 || fabs(bet) > .99*RPT_INF) return bet;

  /* II. check if the distribution is single-sided, which means infinity beta */
  for (hl = hr = 0, i = 0; i < t->m; i++) {
    e = t->emin + i * t->edel;
    if (e < 0) hl += t->h[i];
    else if (e > 0) hr += t->h[i];
  }
  if (hl <= 0.) return RPT_INF;
  if (hr <= 0.) return -RPT_INF;

  /* III. increase |bet|, such that <exp(-bet*de)>  > 1 */
  for (dbet = bet; ; bet *= 1.4)
    if ((f = rpti_getf(t, bet, &df)) > 1) break;
  printf("cnt %d (%d, %d), bet %g, expanded to %g\n", hl+hr, hl, hr, dbet, bet);

  dbmax = fabs(bet*0.1); /* limit the maximal amount of updating */

  /* IV. iteratively refine the temperature */
  for (it = 0; it <= 1000; it++) {
    f = rpti_getf(t, bet, &df);
    /* f should be one, derivative is -df */
    if (fabs(f - 1) < 1e-14 || fabs(df) < 1e-12) break;
    dbet = dblconfine((f - 1)/df, -dbmax, dbmax);
    if (verbose) printf("f %g, df %g, bet %g, dbet %g\n", f, df, bet, dbet);
    bet += dbet;
  }
  if (verbose) printf("\n\n");

  return bet;
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

  xfopen(fp, fn, "w", return -1);
  for (i = i0; i <= i1; i++)
    fprintf(fp, "%d %g %d\n", t->emin + i * t->edel, 1.0*t->h[i]/t->edel, t->h[i]);
  fclose(fp);
  return 0;
}

#endif /* RPT_H__ */

