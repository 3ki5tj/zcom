#include "wl.c"
#ifndef TMH_H__
#define TMH_H__

/* 0: low energy; 1: high energy 
 * if tp0 < tp1: tp is a temperature-like quantity, and rho ~ exp[ -H(erg)/tp ]; 
 * if tp0 > tp1: tp is a beta-like quantity, and rho ~ exp[ -H(erg)*tp ] 
 * */
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
  double scl; /* prefactor in effective Hamiltonian */
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

  wlcvg_t *wl; /* Wang-Landau convergence */
} tmh_t;

tmh_t *tmh_open(double tp0, double tp1, double dtp,
    double erg0, double erg1, double derg,
    double emin, double emax, double de,
    double ensexp, int dhdeorder);
void tmh_close(tmh_t *m);
int tmh_savedhde(tmh_t *m, const char *fn, double amp, double t);
int tmh_loaddhde(tmh_t *m, const char *fn, double *amp, double *t);
int tmh_savetp(tmh_t *m, const char *fn);
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

/* return dH/dE at given energy 'erg', which should be the current potential energy */
INLINE double tmh_getdhde(tmh_t *m, double erg)
{
  double erg0 = m->erg0, erg1 = m->erg1, derg = m->derg, *dhde = m->dhde, x, ix;
  int ergn = m->ergn, ie;

  if (erg <= erg0) return dhde[0];
  else if (erg >= erg1) return dhde[ergn];
  x = modf((erg - erg0) / derg, &ix);
  ie = (int) ix;
  return m->dhdeorder == 0 ? dhde[ie] : dhde[ie] * (1 - x) + dhde[ie + 1] * x;
}

/* update dH/dE curve
 * Note: the location of updating correspond to m->ec instead of erg */
INLINE void tmh_updatedhde(tmh_t *m, double erg, double amp)
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
    TMH_UPDHDE(m->iec + 1, del);
  }
}

/* update dH/dE curve with a diffusive force */
INLINE void tmh_diffusedhde(tmh_t *m, double amp)
{
  die_if (m->dhdeorder == 0, "not for order %d", m->dhdeorder);
  if (m->iec > 0) TMH_UPDHDE(m->iec, amp);
  if (m->iec < m->ergn) TMH_UPDHDE(m->iec + 1, -amp);
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

/* compute dH = H(e1) - H(e0) from integrating the dhde curve */
INLINE double tmh_hdif(tmh_t *m, double e1, double e0)
{
  int ie, iel, ieh, sgn;
  double dh = 0, de, k, el0, eh0, el, eh;

  /* ensure e1 > e0 */
  if (e1 < e0) { eh = e0; el = e1; sgn = -1; }
  else { eh = e1; el = e0; sgn = 1; }

  if (eh < m->erg0) { /* first energy bin */
    dh = (eh - el) * m->dhde[0];
  } else if (el > m->erg1) { /* last energy bin */
    dh = (eh - el) * m->dhde[m->ergn];
  } else {
    /* energy index */
    if (el < m->erg0) {
      dh += (m->erg0 - el) * m->dhde[0];
      el = m->erg0 + 1e-8;
      iel = 0;
    } else {
      iel = (int)((el - m->erg0) / m->derg);
    }
    if (eh >= m->erg1) {
      dh += (eh - m->erg1) * m->dhde[m->ergn];
      eh = m->erg1 - 1e-8;
      ieh = m->ergn - 1;
    } else {
      ieh = (int)((eh - m->erg0) / m->derg);
    }
    if (m->dhdeorder == 0) { /* zeroth order: dhde is constant within a bin */
      if (iel == ieh) {
        dh += (eh - el) * m->dhde[iel];
      } else if (iel < ieh) {
        /* dh at the two terminal energy bin */
        dh += (m->erg0 + (iel+1)*m->derg - el) * m->dhde[iel]
            + (eh - (m->erg0 + m->derg*ieh)) * m->dhde[ieh];
        for (ie = iel+1; ie < ieh; ie++) /* integrate dH/dE */
          dh += m->dhde[ie] * m->derg;
      }
    } else { /* first order: dhde is linear to erg */
      if (iel == ieh) {
        k = (m->dhde[iel+1] - m->dhde[iel]) / m->derg;
        el0 = m->erg0 + iel * m->derg;
        dh += (eh - el) * (m->dhde[iel] + k * (.5f*(el+eh) - el0));
      } else if (iel < ieh) {
        /* dh at the two terminal energy bin */
        el0 = m->erg0 + (iel + 1) * m->derg;
        de = el0 - el;
        k = (m->dhde[iel + 1] - m->dhde[iel]) / m->derg;
        dh += de * (m->dhde[iel + 1] - .5 * k * de);
        eh0 = m->erg0 + m->derg * ieh;
        de = eh - eh0;
        k = (m->dhde[ieh + 1] - m->dhde[ieh]) / m->derg;
        dh += de * (m->dhde[ieh] + .5 * k * de);
        for (ie = iel + 1; ie < ieh; ie++) /* integrate dH/dE */
          dh += .5 * (m->dhde[ie] + m->dhde[ie + 1]) * m->derg;
      }
    }
  }
  return dh * sgn;
}

/* temperature move using a Langevin equation */
#define tmh_lgvmove(m, enow, lgvdt) \
  tmh_langvmove(m, tmh_hdif(m, enow, m->ec), lgvdt)
INLINE int tmh_langvmove(tmh_t *m, double dh, double lgvdt)
{
  double bexp, amp, scl = m->scl;

  if (m->dtp > 0) { /* temperature-like move */
    double tp, tp2, tpl = m->tp0, tph = m->tp1;

    tp = m->tp / scl;
    amp = tp * sqrt(2 * lgvdt);
    bexp = 2. - m->ensexp;
    tp2 = tp + (dh + bexp * tp) * lgvdt + grand0() * amp;
    tp2 *= scl;
    if (tp2 >= tpl && tp2 <= tph) {
      tmh_settp(m, tp2);
      return 1;
    }
  } else { /* beta-like move */
    double bet, bet2, betl = m->tp1, beth = m->tp0;

    bet = m->tp * scl;
    amp = sqrt(2 * lgvdt);
    bet2 = bet - (dh + m->ensexp * bet) * lgvdt + grand0() * amp;
    bet2 /= scl;
    if (bet2 >= betl && bet2 <= beth) {
      tmh_settp(m, bet2);
      return 1;
    }
  }
  return 0;
}

/* initialize amplitude of updating */
INLINE void tmh_initwlcvg(tmh_t *m, double ampmax, double ampfac, double perc,
    double ampc)
{
  double tp0, tp1, dtp;

  if (m->dtp > 0) { tp0 = m->tp0; tp1 = m->tp1; dtp = m->dtp; }
  else { tp0 = m->tp1; tp1 = m->tp0; dtp = -m->dtp; }
  m->wl = wlcvg_open(ampmax, ampfac, perc, ampc, tp0, tp1, dtp);
}

/* easy temperature move
 * famp: tweaking factor */
INLINE int tmh_ezmove(tmh_t *m, double epot, double famp, double lgvdt)
{
  tmh_eadd(m, epot);
  die_if (m->wl == NULL, "call tmh_initamp first, %p\n", (void *) m->wl);
  wlcvg_update(m->wl, m->tp); /* compute updating amplitude */
  tmh_updatedhde(m, epot, m->wl->lnf * famp);
  return tmh_lgvmove(m, epot, lgvdt);
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

INLINE int tmh_save(tmh_t *m, const char *fntp, const char *fnehis,
    const char *fndhde, double amp, double t)
{
  tmh_savetp(m, fntp);
  tmh_savedhde(m, fndhde, amp, t);
  tmh_saveehis(m, fnehis);
  return 0;
}

INLINE int tmh_load(tmh_t *m, const char *fnehis,
    const char *fndhde, double *amp, double *t)
{
  if (tmh_loaddhde(m, fndhde, amp, t) != 0) return -1;
  if (tmh_loadehis(m, fnehis) != 0) return -1;
  return 0;
}

#endif /* TMH_H__ */

