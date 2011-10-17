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
  double *lnz; /* partition function */
  double *lng; /* density of states */
  double *mh; /* modified Hamiltonian  */
} tmh_t;

tmh_t *tmh_open(double tp0, double tp1, double dtp,
    double erg0, double erg1, double derg,
    double emin, double emax, double de,
    double ensexp, int dhdeorder);
void tmh_close(tmh_t *m);
double tmh_hdif(tmh_t *m, double eb, double ea);
int tmh_tlgvmove(tmh_t *m, double enow, double lgvdt);
int tmh_savedhde(tmh_t *m, const char *fn, double amp, double t);
int tmh_loaddhde(tmh_t *m, const char *fn, double *amp, double *t);
int tmh_savetp(tmh_t *m, const char *fn);
int tmh_save(tmh_t *m, const char *fntp, const char *fnehis, 
    const char *fndhde, double amp, double t);
int tmh_load(tmh_t *m, const char *fnehis, 
    const char *fndhde, double *amp, double *t);
int tmh_loaderange(const char *fn, 
    double *tp0, double *tp1, double *dtp,
    double *erg0, double *erg1, double *derg, 
    double *emin, double *emax, double *de,
    double *ensexp, int *dhdeorder);
int tmh_calcdos(tmh_t *m, int itmax, double tol, 
    const char *fndos, const char *fnlnz);

/* set the current temperature */
INLINE void tmh_settp(tmh_t *m, double tp)
{
#ifndef TMH_NOCHECK
  die_if ((m->tp1 - tp) * m->dtp < 0.0 || (tp - m->tp0) * m->dtp < 0.0, 
      "temperature %g not in (%g, %g), dtp %g", tp, m->tp0, m->tp1, m->dtp);
#endif
  m->tp = tp;
  m->itp = (int)((m->tp - m->tp0)/m->dtp);
  m->ec = m->erg0 + (m->tp - m->tp0)*m->dergdt;
  m->iec = (int)((m->ec - m->erg0)/m->derg);
}

/* retrieve local dhde */
INLINE double tmh_getdhde(tmh_t *m, double e, int ie)
{
  if (m->dhdeorder == 0) {
#ifndef TMH_NOCHECK
    die_if (ie < 0 || ie >= m->en, "overflow ie %d en %d\n", ie, m->en);
#endif
    return m->dhde[ie];
  } else {
    double lam = (e - (m->erg0 + ie*m->derg))/m->derg;
#ifndef TMH_NOCHECK
    die_if (lam < 0. || lam > 1., 
        "cannot interpolate, e %g, %d, %g %g %g\n",
        e, ie, m->erg0 + ie*m->derg, m->erg0, m->derg);
#endif
    return m->dhde[ie]*(1-lam) + m->dhde[ie+1]*lam;
  }
}

/* update dhde curve */
INLINE void tmh_dhdeupdate(tmh_t *m, double erg, double amp)
{
  double del = amp * (erg - m->ec);

#ifdef TMH_NOCHECK
  #define TMH_UPDHDE(i, del) m->dhde[i] += del; 
#else
  #define TMH_UPDHDE(i, del) { \
  m->dhde[i] += del; \
  if (m->dhde[i] < m->dhdemin) \
    m->dhde[i] = m->dhdemin; \
  else if (m->dhde[i] > m->dhdemax) \
    m->dhde[i] = m->dhdemax; }
#endif

  if (m->dhdeorder == 0) {
    TMH_UPDHDE(m->iec, del);
    if (m->iec == m->ergn - 1) /* last bin */
      m->dhde[m->ergn] = m->dhde[m->iec];
  } else {
    del *= .5;
    TMH_UPDHDE(m->iec, del);
    TMH_UPDHDE(m->iec+1, del);
  }
}

INLINE void tmh_eadd(tmh_t *m, double erg)
{
  int ie;
#ifndef TMH_NOCHECK
  if (erg < m->emin || erg > m->emax) return;
#endif
  ie = (int)((erg - m->emin)/m->de);
#ifndef TMH_NOCHECK
  die_if (ie < 0 || ie >= m->en, "ie = %d, erg %g emin %g de %g output range\n", 
      ie, erg, m->emin, m->de);  
  die_if (m->itp > m->tpn, "itp = %d, tpn = %d, tp = %g, dtp = %g\n", 
      m->itp, m->tpn, m->tp, m->dtp);
#endif
  m->tpehis[m->itp*m->en + ie] += 1.;
}

INLINE int tmh_saveehis(tmh_t *m, const char *fn)
{
  return histsave(m->tpehis, m->tpn, 
      m->en, m->emin, m->de, HIST_ADDAHALF|HIST_OVERALL, fn);
}

INLINE int tmh_loadehis(tmh_t *m, const char *fn)
{
  return histload(m->tpehis, m->tpn, 
      m->en, m->emin, m->de, HIST_ADDAHALF, fn);
}


#endif /* TMH_H__ */

