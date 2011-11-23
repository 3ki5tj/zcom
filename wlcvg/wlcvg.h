#include "util.h"
#include "hist.c"
#ifndef WLCVG_H__
#define WLCVG_H__
#include <stdio.h>

typedef struct {
  double lnf0;  /* lnf of the first stage, also the maximal value */
  double lnffac; /* factor to multiply on passing to another stage */
  double lnfc;  /* lnf to be used in formula */
  double lnfmin; /* minimal lnf */
  double lnfwl; /* lnf as given by Wang-Landau algorithm */
  double lnfinvt; /* lnf as given value given by lnfc/nupd; */
  double lnf; /* lnf adjusted by both wl and invt */
  double perc, percutoff;
  double nupd; /* the total number of updates */
  int stage; /* stage id */
  int nstcheck; /* frequency of checking histogram flatness */
  unsigned flags;
  hist_t *hs;
} wlcvg_t;

#define WLCVG_UPDLNFC 0x0010

/* to use a pure Wang-Landau scheme, set lnfc = 0
 * to use a pure formula scheme, set percutoff = 0 */
INLINE wlcvg_t *wlcvg_open(double lnfc, double lnf0, double lnffac, double percutoff, double lnfmin,
    double xmin, double xmax, double dx)
{
  wlcvg_t *wl;

  xnew(wl, 1);
  wl->hs = hs_open(1, xmin, xmax, dx);
  wl->lnf = wl->lnfwl = wl->lnf0 = lnf0;
  wl->lnffac = lnffac;
  wl->lnfc = lnfc;
  wl->lnfmin = lnfmin;
  wl->percutoff = percutoff;
  wl->perc = 1.0;
  wl->nupd = 0.0;
  wl->stage = 0;
  wl->nstcheck = 10 * ((int)((xmax - xmin)/dx + 0.99999999));
  wl->flags = 0u;
  return wl;
}

INLINE void wlcvg_close(wlcvg_t *wl)
{
  hs_close(wl->hs);
  free(wl);
}

/* compute histogram flatness */
INLINE double wlcvg_hsflatness(wlcvg_t *wl)
{
  double hmin = 1.0e30, hmax = 0.0, h;
  int i;

  for (i = 0; i < wl->hs->n; i++) {
    h = wl->hs->arr[i];
    if (h < hmin) hmin = h;
    else if (h > hmax) hmax = h;
  }
  return wl->perc = (hmax > 0 ? (hmax - hmin)/(hmax + hmin) : 1.0);
}

/* register x to the histogram */
INLINE void wlcvg_add(wlcvg_t *wl, double x)
{
  wl->nupd += 1.0;
  hs_add(wl->hs, &x, 1, 0);
}

INLINE void wlcvg_shiftstage(wlcvg_t *wl)
{
  wl->stage++;
  wl->perc = 1.0;
  if ((wl->flags & WLCVG_UPDLNFC) && wl->lnfinvt < wl->lnfwl)
    wl->lnfc = wl->nupd * wl->lnfwl;
  wl->lnfwl *= wl->lnffac;
  hs_clear(wl->hs);
}

/* easy driver for updating lnf, return lnf */
INLINE double wlcvg_update(wlcvg_t *wl, double x)
{
  wlcvg_add(wl, x);
  if (fmod(wl->nupd, wl->nstcheck) < 0.1) /* compute histogram flatness */
    wlcvg_hsflatness(wl);
  wl->lnfinvt = wl->lnfc / wl->nupd; /* value given by formula */
  if (wl->percutoff > 0 && wl->lnfinvt < wl->lnfwl) { /* consult the Wang-Landau scheme */
    wl->lnf = wl->lnfwl;
    if (wl->perc < wl->percutoff) wlcvg_shiftstage(wl);
  } else {
    wl->lnf = dblmin(wl->lnfinvt, wl->lnf0);
  }
  return wl->lnf = dblmax(wl->lnf, wl->lnfmin);
}

INLINE void wlcvg_setnupd(wlcvg_t *wl, double nupd, double lnfwl, double lnfc)
{
  modf(nupd, &wl->nupd); /* truncate to the nearest integer */
  if (lnfwl > 0) wl->lnfwl = lnfwl;
  if (lnfc > 0) wl->lnfc = lnfc;
}
#endif /* WLCVG_H__ */
