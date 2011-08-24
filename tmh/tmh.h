#ifndef TMH_H__
#define TMH_H__

typedef struct {
  double tp, ec; /* current temperature, and the expected energy there */
  int itp, iec; /* indices of tp and ec */
  double tp0, tp1, dtp; /* temperature range */
  int tpn; /* number of temperature */
  double emin, emax; /* energy range */
  double de; /* bin size of energy histogram */
  int en; /* number of energy bins */
  double erg0, erg1; /* energy range (erg0, erg1) */
  double derg; /* bin size for the updating energy range */
  int ergn; /* number of the updating energy bins */
  double dergdt; /* (erg1 - erg0)/(tp1 - tp0) */
  int dhdeorder; /* order of dhde interpolation */
  double dhdemin; /* minimal of dhde */
  double dhdemax; /* maximal of dhde */
  double *dhde; /* dH / dE - 1 */
  double *tpehis; /* multipl-temperature energy histogram */
  double ensexp; /* w(T) = 1/T^ensexp */
} tmh_t;

tmh_t *tmh_open(double tp0, double tp1,
    double erg0, double erg1, double derg,
    double emin, double emax, double de,
    double ensexp, int dhdeorder);
void tmh_close(tmh_t *tmh);
double tmh_hdif(tmh_t *tmh, double eb, double ea);
int tmh_tlgvmove(tmh_t *tmh, double enow, double lgvdt);
int tmh_savedhde(tmh_t *tmh, const char *fn, double amp, double t);
int tmh_loaddhde(tmh_t *tmh, const char *fn, double *amp, double *t);
int tmh_savetp(tmh_t *tmh, const char *fn);
int tmh_save(tmh_t *tmh, const char *fntp, const char *fnehis, 
    const char *fndhde, double amp, double t);
int tmh_load(tmh_t *tmh, const char *fnehis, 
    const char *fndhde, double *amp, double *t);
int tmh_loaderange(const char *fn, double *erg0, double *erg1, double *derg, 
    double *emin, double *emax, double *de);

/* set the current temperature */
ZCINLINE void tmh_settp(tmh_t *tmh, double tp)
{
#ifndef TMH_NOCHECK
  die_if (tp > tmh->tp1 || tp < tmh->tp0, "temperature %g not in(%g, %g)", tp, tmh->tp0, tmh->tp1);
#endif
  tmh->tp = tp;
  tmh->itp = (int)((tmh->tp - tmh->tp0)/tmh->dtp);
  tmh->ec = tmh->erg0 + (tmh->tp - tmh->tp0)*tmh->dergdt;
#ifndef TMH_NOCHECK
  tmh->iec = (int)((tmh->ec - tmh->erg0)/tmh->derg);
#else
  tmh->iec = tmh->itp;
#endif
}

/* retrieve local dhde */
ZCINLINE double tmh_getdhde(tmh_t *tmh)
{
  if (tmh->dhdeorder == 0) {
    return tmh->dhde[tmh->iec];
  } else {
    double lam = (tmh->ec - (tmh->erg0 + tmh->iec*tmh->derg))/tmh->derg;
    return tmh->dhde[tmh->iec]*(1-lam) + tmh->dhde[tmh->iec+1]*lam;
  }
}

/* update dhde curve */
ZCINLINE void tmh_dhdeupdate(tmh_t *tmh, double erg, double amp)
{
  double del = amp * (erg - tmh->ec);

#ifdef TMH_NOCHECK
  #define TMH_UPDHDE(i, del) \
    tmh->dhde[i] += del
#else
  #define TMH_UPDHDE(i, del) \
  tmh->dhde[i] += del; \
  if (tmh->dhde[i] < tmh->dhdemin) \
    tmh->dhde[i] = tmh->dhdemin; \
  else if (tmh->dhde[i] > tmh->dhdemax) \
    tmh->dhde[i] = tmh->dhdemax;
#endif

  if (tmh->dhdeorder == 0) {
    TMH_UPDHDE(tmh->iec, del);
    if (tmh->iec == tmh->ergn - 1) /* last bin */
      tmh->dhde[tmh->ergn] = tmh->dhde[tmh->iec];
  } else {
    del *= .5;
    TMH_UPDHDE(tmh->iec, del);
    TMH_UPDHDE(tmh->iec+1, del);
  }
}

ZCINLINE void tmh_eadd(tmh_t *tmh, double erg)
{
  int ie;
#ifndef TMH_NOCHECK
  if (erg < tmh->emin || erg > tmh->emax) return;
#endif
  ie = (int)((erg - tmh->emin)/tmh->de);
#ifndef TMH_NOCHECK
  die_if (ie < 0 || ie >= tmh->en, "ie = %d, erg %g emin %g de %g output range\n", 
      ie, erg, tmh->emin, tmh->de);  
  die_if (tmh->itp > tmh->tpn, "itp = %d, tpn = %d, tp = %g, dtp = %g\n", 
      tmh->itp, tmh->tpn, tmh->tp, tmh->dtp);
#endif
  tmh->tpehis[tmh->itp*tmh->en + ie] += 1.;
}

ZCINLINE int tmh_saveehis(tmh_t *tmh, const char *fn)
{
  return histsave(tmh->tpehis, tmh->tpn, 
      tmh->en, tmh->emin, tmh->de, HIST_ADDAHALF|HIST_OVERALL, fn);
}

ZCINLINE int tmh_loadehis(tmh_t *tmh, const char *fn)
{
  return histload(tmh->tpehis, tmh->tpn, 
      tmh->en, tmh->emin, tmh->de, HIST_ADDAHALF, fn);
}


#endif /* TMH_H__ */

