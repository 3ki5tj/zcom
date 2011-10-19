#include "util.h"
#include "hist.c"
#include "rng.c"

#ifndef TMH_C__
#define TMH_C__

/* tempering with modified Hamiltonian */
#include "tmh.h"

tmh_t *tmh_open(double tp0, double tp1, double dtp,
    double erg0, double erg1, double derg,
    double emin, double emax, double de,
    double ensexp, int dhdeorder)
{
  tmh_t *m;
  int i;

  xnew(m, 1);

  m->scl = 1; /* Hamiltonian scaling factor */

  /* energy histogram range */
  m->de = de;
  die_if(emin >= emax, "Error: emin %g >= emax %g\n", emin, emax);
  m->emin = emin;
  m->emax = emin + dblround(emax - emin, de);
  m->en = (int)((m->emax - m->emin)/de + .5); /* number of energy bins */

  /* the updating energy range */
  die_if(erg0 >= erg1, "Error: erg0 %g >= erg1 %g\n", erg0, erg1);
  m->derg = derg;
  m->erg0 = erg0;
  m->erg1 = erg0 + dblround(erg1 - erg0, derg);
  m->ergn = (int)((m->erg1 - m->erg0)/m->derg + .5);

  /* dhde parameters */
  m->dhdeorder = dhdeorder;
  m->dhdemin = 0.1;
  m->dhdemax = 10.0;
  xnew(m->dhde, m->ergn + 1);
  for (i = 0; i <= m->ergn; i++)
    m->dhde[i] = 1.;

  die_if((tp1 - tp0)*dtp < 0, "Error: tp0 %g, tp1 %g, dtp %g\n", tp0, tp1, dtp);
  if (fabs(dtp) > 0) { /* dtp is explicitly specified */
    m->dtp = dtp;
    m->tp0 = dblround(tp0, dtp); /* round tp0, tp1 to multiples of dtp */
    m->tp1 = dblround(tp1, dtp);
    m->tpn = (int)((m->tp1 - m->tp0)/m->dtp + .5);
  } else {
    m->tp0 = tp0;
    m->tp1 = tp1;
    m->tpn = m->ergn;
    m->dtp = (m->tp1 - m->tp0)/m->tpn * (1. + 1e-12);
  }
  xnew(m->tpehis, m->tpn*m->en);

  m->ensexp = ensexp;
  tmh_settp(m, (m->tp0+m->tp1)*.5);
  m->dergdt = (m->erg1 - m->erg0)/(m->tp1 - m->tp0);

  xnew(m->lnz, m->tpn);
  xnew(m->lng, m->en);
  xnew(m->mh, m->en+1);

  return m;
}

void tmh_close(tmh_t *m)
{
  if (m != NULL) {
    free(m->dhde);
    free(m->tpehis);
    free(m->lnz);
    free(m->lng);
    free(m->mh);
    free(m);
  }
}

/* write dhde and overall energy distribution */
int tmh_savedhde(tmh_t *m, const char *fn, double amp, double t)
{
  int ie;
  FILE *fp;

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# 2 %g %g %d %g %g %d %g %g %d %g %d %g %g %g\n",
      m->erg0, m->derg, m->ergn,
      m->emin, m->de,   m->en,
      m->tp0,  m->dtp,  m->tpn,
      m->ensexp, m->dhdeorder,
      amp, t, m->tp);
  for (ie = 0; ie <= m->ergn; ie++) {
    fprintf(fp, "%g %g\n", m->erg0 + ie*m->derg, m->dhde[ie]);
  }
  fclose(fp);
  return 0;
}

/* read dhde and overall energy distribution */
int tmh_loaddhde(tmh_t *m, const char *fn, double *amp, double *t)
{
  int ie, ergn, en, ver, tpn, dhdeorder, next;
  FILE *fp;
  char s[1024], *p;
  double emin, de, derg, erg, erg0, tp0, dtp, dhde, ensexp;

  xfopen(fp, fn, "r", return -1);
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "cannot read the first line %s\n", fn);
    goto ERR;
  }
  if (7 != sscanf(s, " # %d%lf%lf%d%lf%lf%d%n",
        &ver, &erg0, &derg, &ergn, &emin, &de, &en, &next)) {
    fprintf(stderr, "corrupted info line %s", s);
    goto ERR;
  }
  if (ver < 1 ||
      ergn != m->ergn || fabs(derg - m->derg) > 1e-5 || fabs(erg0 - m->erg0) > 1e-5 ||
      en != m->en || fabs(de - m->de) > 1e-5 || fabs(emin - m->emin) > 1e-5) {
    fprintf(stderr, "bad energy ergn %d %d, erg0 %g %g, derg %g %g, en %d %d, emin %g %g, de %g %g\n",
        ergn, m->ergn, erg0, m->erg0, derg, m->derg,
        en, m->en, emin, m->emin, de, m->de);
    goto ERR;
  }
  p = s + next;
  if (ver >= 2) {
    if (5 != sscanf(p, "%lf%lf%d%lf%d%n", &tp0, &dtp, &tpn, &ensexp, &dhdeorder, &next)) {
      fprintf(stderr, "corrupted info line p2 %s", s);
      goto ERR;
    }
    p += next;
    if (tpn != m->tpn || fabs(tp0 - m->tp0) > 1e-6 || fabs(dtp - m->dtp) > 1e-6) {
      fprintf(stderr, "bad temperature range tpn %d %d, tp0 %g %g, dtp %g %g\n",
          tpn, m->tpn, tp0, m->tp0, dtp, m->dtp);
      goto ERR;
    }
  }
  if (3 != sscanf(p, "%lf%lf%lf", amp, t, &m->tp)) {
    fprintf(stderr, "corrupted info line p3 %s", s);
    goto ERR;
  }

  for (ie = 0; ie < m->ergn; ie++) {
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "cannot read line %d\n", ie);
      goto ERR;
    }
    if (2 != sscanf(s, "%lf%lf", &erg, &dhde)) {
      fprintf(stderr, "cannot read energy and dhde at %d\n, s = %s", ie, s);
      goto ERR;
    }
    erg0 = m->erg0 + ie*m->derg;
    if (fabs(erg0 - erg) > m->derg*.1) {
      fprintf(stderr, "energy %g, should be %g\n", erg, erg0);
      goto ERR;
    }
    m->dhde[ie] = dhde;
  }
  fclose(fp);
  return 0;
ERR:
  fclose(fp);
  return -1;
}

/* get energy range from dhde file */
int tmh_loaderange(const char *fn,
    double *tp0, double *tp1, double *dtp,
    double *erg0, double *erg1, double *derg,
    double *emin, double *emax, double *de,
    double *ensexp, int *dhdeorder)
{
  int en, ergn, tpn, ver, next;
  FILE *fp;
  char s[1024];

  xfopen(fp, fn, "r", return -1);
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "cannot read the first line %s\n", fn);
    goto ERR;
  }

  if (7 != sscanf(s, " # %d%lf%lf%d%lf%lf%d%n",
        &ver, erg0, derg, &ergn, emin, de, &en, &next) || ver < 1) {
    fprintf(stderr, "corrupted info line %s", s);
    goto ERR;
  }
  *erg1 = *erg0 + ergn*(*derg);
  *emax = *emin + en*(*de);

  if (ver >= 2) { /* additional information */
    if (5 != sscanf(s+next, "%lf%lf%d%lf%d", tp0, dtp, &tpn, ensexp, dhdeorder)) {
      fprintf(stderr, "corrupted info line %s", s);
      goto ERR;
    }
    *tp1 = *tp0 + tpn * (*dtp);
  }

  fclose(fp);
  return 0;
ERR:
  fclose(fp);
  return -1;
}

/* write temperature histogram */
int tmh_savetp(tmh_t *m, const char *fn)
{
  int i, j;
  double *eh, erg, cnt, esm, e2sm, eav, edv;
  FILE *fp;

  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %g %g %d\n", m->tp0, m->dtp, m->tpn);
  for (i = 0; i < m->tpn; i++) {
    eh = m->tpehis + i*m->en;
    for (cnt = esm = e2sm = 0., j = 0; j < m->en; j++) {
      erg = m->emin + (j + .5) * m->de;
      cnt += eh[j];
      esm += eh[j]*erg;
      e2sm += eh[j]*(erg*erg);
    }
    if (cnt > 1e-8) {
      eav = esm/cnt;
      edv = sqrt(e2sm/cnt - eav*eav);
    } else {
      eav = edv = 0.;
    }
    fprintf(fp, "%g %g %g %g\n",
        m->tp0 + (i+.5)*m->dtp, cnt, eav, edv);
  }
  fclose(fp);
  return 0;
}

/* calculate the modified Hamiltonian */
static int tmh_calcmh(tmh_t *m)
{
  int i, en = m->en;
  double erg, hm, dh, de = m->de, emin = m->emin;

  hm = emin;
  for (i = 0; i < en; i++) {
    erg = emin + (i + .5) * de;
    dh = de * tmh_getdhde(m, erg);
    m->mh[i] = hm + .5 * dh;
    hm += dh;
  }
  return 0;
}

/* iteratively compute the density of states */
int tmh_calcdos(tmh_t *m, int itmax, double tol,
    const char *fndos, const char *fnlnz)
{
  int i, j, it, ie0, ie1, en = m->en, tpn = m->tpn;
  double x, dif, lnz0, db;
  double *lnn, *lnm, *bet, *lnz1;
  const double LOG0 = -1e10;
  FILE *fp;

  if (itmax <= 0) itmax = 1000;

  /* determine nonempty energy range for more efficient loops */
  for (ie0 = 0; ie0 < en; ie0++) {
    for (j = 0; j < m->tpn; j++)
      if (m->tpehis[j*en + ie0] > 0.) break;
    if (j < m->tpn) break;
  }
  for (ie1 = en; ie1 > ie0; ie1--) {
    for (j = 0; j < m->tpn; j++)
      if (m->tpehis[j*en + ie1-1] > 0.) break;
    if (j < m->tpn) break;
  }

  /* n[j] is the total number of visits to temperature j
   * m[i] is the total number of visits to energy i */
  xnew(lnn, tpn);
  xnew(lnm, en);
  for (i = 0; i < en; i++) lnm[i] = 0.;
  for (j = 0; j < tpn; j++) {
    lnn[j] = 0.;
    for (i = ie0; i < ie1; i++) {
      x = m->tpehis[j*en + i];
      lnn[j] += x;
      lnm[i] += x;
    }
  }
  for (j = 0; j < tpn; j++)
    lnn[j] = (lnn[j] > 0.) ? log(lnn[j]) : LOG0;
  for (i = 0; i < en; i++)
    lnm[i] = (lnm[i] > 0.) ? log(lnm[i]) : LOG0;

  xnew(bet, tpn);
  for (j = 0; j < tpn; j++) {
    double tp = m->tp0 + (j + .5) * m->dtp;
    if (m->dtp > 0.) bet[j] = m->scl / tp;
    else bet[j] = m->scl * tp;
  }

  /* get mh and lnz */
  for (i = 0; i < en; i++) m->lng[i] = LOG0;
  tmh_calcmh(m);
  /* estimate initial lnz */
  for (lnz0 = 0., j = 0; j < tpn; j++) {
    x = m->erg0 + (j+.5)*m->dtp * m->dergdt;
    i = (int)((x - m->emin)/m->de);
    die_if (i < 0 || i >= en, "i %d, x %g\n", i, x);
    db = 1.0/(m->tp0 + (j+1)*m->dtp)
       - 1.0/(m->tp0 + j*m->dtp);
    x = m->mh[i]*db;
    m->lnz[j] = lnz0 - x*.5;
    lnz0 -= x;
  }
  for (j = tpn - 1; j >= 0; j--)
     m->lnz[j] -= m->lnz[0];

  xnew(lnz1, tpn);

  /* repeat until convergence */
  for (it = 0; it < itmax; it++) {
    /* compute the density of states */
    for (i = ie0; i < ie1; i++) {
      for (x = LOG0, j = 0; j < tpn; j++) {
        x = lnadd(x, lnn[j] - bet[j]*m->mh[i] - m->lnz[j]);
      }
      m->lng[i] = (lnm[i] < LOG0+.1) ? LOG0 : (lnm[i] - x);
    }

    /* update partition function */
    for (j = 0; j < tpn; j++) {
      for (x = LOG0, i = ie0; i < ie1; i++)
        x = lnadd(x, m->lng[i] - bet[j]*m->mh[i]);
      lnz1[j] = x;
    }
    for (j = tpn - 1; j >= 0; j--)
      lnz1[j] -= lnz1[0];

    /* check difference */
    for (dif = 0., j = 1; j < tpn; j++) {
      x = fabs(m->lnz[j] - lnz1[j]);
      if (x > dif) dif = x;
    }
    for (j = 0; j < tpn; j++) m->lnz[j] = lnz1[j];
    if (dif < tol) break;
  }

  /* write dos */
  if (fndos && (fp = fopen(fndos, "w")) != NULL) {
    for (i = ie0; i < ie1; i++) {
      x = m->emin + (i + .5)*m->de;
      fprintf(fp, "%g %g %g %g\n", x, m->lng[i] - m->lng[ie0],
          exp(lnm[i]), m->mh[i]);
    }
    fclose(fp);
  }

  /* write lnz */
  if (fnlnz && (fp = fopen(fnlnz, "w")) != NULL) {
    for (j = 0; j < tpn; j++) {
      fprintf(fp, "%g %g %g\n", bet[j], m->lnz[j], exp(lnn[j]));
    }
    fclose(fp);
  }
  free(lnm);
  free(lnn);
  free(bet);
  free(lnz1);
  return 0;
}

#endif

