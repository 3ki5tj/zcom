#include "util.c"
#include "hist.c"
#include "rng.c"

#ifndef TMH_C__
#define TMH_C__

/* tempering with modified Hamiltonian */
#include "tmh.h"

/* 0: cold; 1: hot */
tmh_t *tmh_open(double tp0, double tp1, double dtp, 
    double erg0, double erg1, double derg,
    double emin, double emax, double de,
    double ensexp, int dhdeorder)
{
  tmh_t *tmh;
  int i;

  xnew(tmh, 1);

  /* energy histogram range */
  tmh->de = de;
  die_if(emin >= emax, "Error: emin %g >= emax %g\n", emin, emax);
  tmh->emin = emin;
  tmh->emax = emin + dblround(emax - emin, de);
  tmh->en = (int)((tmh->emax - tmh->emin)/de + .5); /* number of energy bins */

  /* the updating energy range */
  die_if(erg0 >= erg1, "Error: erg0 %g >= erg1 %g\n", erg0, erg1);
  tmh->derg = derg;
  tmh->erg0 = erg0;
  tmh->erg1 = erg0 + dblround(erg1 - erg0, derg);
  tmh->ergn = (int)((tmh->erg1 - tmh->erg0)/tmh->derg + .5);

  /* dhde parameters */ 
  tmh->dhdeorder = dhdeorder; 
  tmh->dhdemin = 0.1;
  tmh->dhdemax = 10.0;
  xnew(tmh->dhde, tmh->ergn + 1);
  for (i = 0; i <= tmh->ergn; i++)
    tmh->dhde[i] = 1.;

  die_if(tp0 >= tp1, "Error: T0 %g >= T1 %g\n", tp0, tp1);
  if (dtp <= 0) {
    tmh->tp0 = tp0;
    tmh->tp1 = tp1;
    tmh->tpn = tmh->ergn;
    tmh->dtp = (tmh->tp1 - tmh->tp0)/tmh->tpn * (1. + 1e-12);
  } else {
    tmh->dtp = dtp;
    tmh->tp0 = dblround(tp0, dtp);
    tmh->tp1 = dblround(tp1, dtp);
    tmh->tpn = (int)((tmh->tp1 - tmh->tp0)/tmh->dtp + .5);
  }
  xnew(tmh->tpehis, tmh->tpn*tmh->en);

  tmh->ensexp = ensexp;
  tmh_settp(tmh, (tmh->tp0+tmh->tp1)*.5);
  tmh->dergdt = (tmh->erg1 - tmh->erg0)/(tmh->tp1 - tmh->tp0);

  xnew(tmh->lnz, tmh->tpn);
  xnew(tmh->lng, tmh->en);
  xnew(tmh->mh, tmh->en+1);

  return tmh;
}

void tmh_close(tmh_t *tmh)
{
  if (tmh != NULL) {
    free(tmh->dhde);
    free(tmh->tpehis);
    free(tmh->lnz);
    free(tmh->lng);
    free(tmh->mh);
    free(tmh);
  }
}

static double tmh_hdif0(tmh_t *tmh, double eh, double el)
{
  int ie, iel, ieh;
  double dh;

  if (eh < tmh->erg0) { /* first energy bin */
    dh = (eh - el)*tmh->dhde[0];
  } else if (el > tmh->erg1) { /* last energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ergn];
  } else {
    dh = 0.;
    /* overhead index */
    if (el < tmh->erg0) {
      dh += (tmh->erg0 - el)*tmh->dhde[0];
      el = tmh->erg0 + 1e-8;
    }
    if (eh >= tmh->erg1) {
      dh += (eh - tmh->erg1)*tmh->dhde[tmh->ergn];
      eh = tmh->erg1 - 1e-8;
    }
    /* energy index */
    iel = (int)((el - tmh->erg0)/tmh->derg);
    ieh = (int)((eh - tmh->erg0)/tmh->derg);
    if (iel == ieh) {
      dh += (eh - el)*tmh->dhde[iel];
    } else if (iel < ieh) {
      /* dh at the two terminal energy bin */
      dh += (tmh->erg0 + (iel+1)*tmh->derg - el) * tmh->dhde[iel]
          + (eh - (tmh->erg0 + tmh->derg*ieh)) * tmh->dhde[ieh];
      /* integrate dH/dE */
      for (ie = iel+1; ie < ieh; ie++)
        dh += tmh->dhde[ie]*tmh->derg;
    }
  }
  return dh; 
}

static double tmh_hdif1(tmh_t *tmh, double eh, double el)
{
  int ie, iel, ieh;
  double dh, de, k, el0, eh0;

  if (eh < tmh->erg0) { /* first energy bin */
    dh = (eh - el)*tmh->dhde[0];
  } else if (el > tmh->erg1) { /* last energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ergn];
  } else {
    dh = 0.;
    if (el < tmh->erg0) {
      dh += (tmh->erg0 - el)*tmh->dhde[0];
      el = tmh->erg0 + 1e-8;
    }
    if (eh >= tmh->erg1) {
      dh += (eh - tmh->erg1)*tmh->dhde[tmh->ergn];
      eh = tmh->erg1 - 1e-8;
    }
    /* energy index */
    iel = (int)((el - tmh->erg0)/tmh->derg);
    ieh = (int)((eh - tmh->erg0)/tmh->derg);
    if (iel == ieh) {
      k = (tmh->dhde[iel+1] - tmh->dhde[iel])/tmh->derg;
      el0 = tmh->erg0 + iel*tmh->derg;
      dh += (eh - el)*(tmh->dhde[iel] + k * (.5f*(el+eh) - el0));
    } else if (iel < ieh) {
      /* dh at the two terminal energy bin */
      el0 = tmh->erg0 + (iel+1)*tmh->derg; 
      de = el0 - el;
      k = (tmh->dhde[iel+1] - tmh->dhde[iel])/tmh->derg;
      dh += de * (tmh->dhde[iel+1] - k*.5*de);
      eh0 = tmh->erg0 + tmh->derg*ieh;
      de = eh - eh0;
      k = (tmh->dhde[ieh+1] - tmh->dhde[ieh])/tmh->derg;
      dh += de * (tmh->dhde[ieh] + k*.5*de);
      /* integrate dH/dE */
      for (ie = iel+1; ie < ieh; ie++)
        dh += .5*(tmh->dhde[ie] + tmh->dhde[ie+1])*tmh->derg;
    }
  }
  return dh; 
}

/* d(H) = H(e1) - H(e0), */
double tmh_hdif(tmh_t *tmh, double e1, double e0)
{
  int sgn;
  double tmp;

  /* to make sure e1 > e0 */
  if (e1 < e0) {
    sgn = -1;
    tmp = e1, e1 = e0, e0 = tmp;
  } else sgn = 1;

  return sgn * ((tmh->dhdeorder == 0) ?
    tmh_hdif0(tmh, e1, e0) : tmh_hdif1(tmh, e1, e0));
}

/* temperature move using a Langevin equation */
int tmh_tlgvmove(tmh_t *tmh, double enow, double lgvdt)
{
  double derg, tpp, bexp = 2.-tmh->ensexp, amp;
  int lgvtype = 1;
  derg = tmh_hdif(tmh, enow, tmh->ec);
  if (lgvtype == 0) {
    amp = sqrt(2*lgvdt);
    tpp = 1.0 / ( 1.0/(tmh->tp) - (derg+bexp*tmh->tp)*lgvdt + grand0()*amp );
  } else {
    amp = tmh->tp*sqrt(2*lgvdt);
    tpp = tmh->tp + (derg+bexp*tmh->tp)*lgvdt + grand0()*amp;
  }
  if (tpp > tmh->tp0 && tpp < tmh->tp1) {
    tmh_settp(tmh, tpp);
    return 1;
  }
  return 0;
}

/* write dhde and overall energy distribution */
int tmh_savedhde(tmh_t *tmh, const char *fn, double amp, double t)
{
  int ie; 
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
  fprintf(fp, "# 2 %g %g %d %g %g %d %g %g %d %g %d %g %g %g\n", 
      tmh->erg0, tmh->derg, tmh->ergn,
      tmh->emin, tmh->de,   tmh->en,
      tmh->tp0,  tmh->dtp,  tmh->tpn,
      tmh->ensexp, tmh->dhdeorder,
      amp, t, tmh->tp);
  for (ie = 0; ie <= tmh->ergn; ie++) {
    fprintf(fp, "%g %g\n", tmh->erg0 + ie*tmh->derg, tmh->dhde[ie]);
  }
  fclose(fp);
  return 0;
}

/* read dhde and overall energy distribution */
int tmh_loaddhde(tmh_t *tmh, const char *fn, double *amp, double *t)
{
  int ie, ergn, en, ver, tpn, dhdeorder, next;
  FILE *fp;
  char s[1024], *p;
  double emin, de, derg, erg, erg0, tp0, dtp, dhde, ensexp;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
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
      ergn != tmh->ergn || fabs(derg - tmh->derg) > 1e-5 || fabs(erg0 - tmh->erg0) > 1e-5 ||
      en != tmh->en || fabs(de - tmh->de) > 1e-5 || fabs(emin - tmh->emin) > 1e-5) {
    fprintf(stderr, "bad energy ergn %d %d, erg0 %g %g, derg %g %g, en %d %d, emin %g %g, de %g %g\n",
        ergn, tmh->ergn, erg0, tmh->erg0, derg, tmh->derg,
        en, tmh->en, emin, tmh->emin, de, tmh->de);
    goto ERR;
  }
  p = s + next;
  if (ver >= 2) {
    if (5 != sscanf(p, "%lf%lf%d%lf%d%n", &tp0, &dtp, &tpn, &ensexp, &dhdeorder, &next)) {
      fprintf(stderr, "corrupted info line p2 %s", s);
      goto ERR;
    }
    p += next;
    if (tpn != tmh->tpn || fabs(tp0 - tmh->tp0) > 1e-6 || fabs(dtp - tmh->dtp) > 1e-6) {
      fprintf(stderr, "bad temperature range tpn %d %d, tp0 %g %g, dtp %g %g\n",
          tpn, tmh->tpn, tp0, tmh->tp0, dtp, tmh->dtp);
      goto ERR;
    }
  }
  if (3 != sscanf(p, "%lf%lf%lf", amp, t, &tmh->tp)) {
    fprintf(stderr, "corrupted info line p3 %s", s);
    goto ERR;
  }

  for (ie = 0; ie < tmh->ergn; ie++) {
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "cannot read line %d\n", ie);
      goto ERR;
    }
    if (2 != sscanf(s, "%lf%lf", &erg, &dhde)) {
      fprintf(stderr, "cannot read energy and dhde at %d\n, s = %s", ie, s);
      goto ERR;
    }
    erg0 = tmh->erg0 + ie*tmh->derg;
    if (fabs(erg0 - erg) > tmh->derg*.1) {
      fprintf(stderr, "energy %g, should be %g\n", erg, erg0);
      goto ERR;
    }
    tmh->dhde[ie] = dhde;
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

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
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
int tmh_savetp(tmh_t *tmh, const char *fn)
{
  int i, j;
  double *eh, erg, cnt, esm, e2sm, eav, edv;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %g %g %d\n", tmh->tp0, tmh->dtp, tmh->tpn);
  for (i = 0; i < tmh->tpn; i++) {
    eh = tmh->tpehis + i*tmh->en;
    for (cnt = esm = e2sm = 0., j = 0; j < tmh->en; j++) {
      erg = tmh->emin + (j + .5) * tmh->de;
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
        tmh->tp0 + (i+.5)*tmh->dtp, cnt, eav, edv);
  }
  fclose(fp);
  return 0;
}

int tmh_save(tmh_t *tmh, const char *fntp, const char *fnehis, 
    const char *fndhde, double amp, double t)
{
  tmh_savetp(tmh, fntp);
  tmh_savedhde(tmh, fndhde, amp, t);
  tmh_saveehis(tmh, fnehis);
  return 0;
}

int tmh_load(tmh_t *tmh, const char *fnehis, 
    const char *fndhde, double *amp, double *t)
{
  if (tmh_loaddhde(tmh, fndhde, amp, t) != 0) return -1;
  if (tmh_loadehis(tmh, fnehis) != 0) return -1;
  return 0;
}

/* calculate the modified Hamiltonian */
static int tmh_calcmh(tmh_t *tmh)
{
  int i, ie;
  double erg, hm, dh;

  hm = tmh->emin;
  for (i = 0; i < tmh->en; i++) {
    erg = tmh->emin + (i+.5)*tmh->de;
    if (erg <= tmh->erg0) {
      dh = tmh->dhde[0];
    } else if (erg >= tmh->erg1) {
      dh = tmh->dhde[tmh->ergn];
    } else {
      ie = (int)((erg - tmh->erg0)/tmh->derg);
      die_if (ie < 0 || ie >= tmh->ergn,
          "ie %d, erg %g, erg0 %g, derg %g",
          ie, erg, tmh->erg0, tmh->derg);
      dh = tmh_getdhde(tmh, erg, ie);
    }
    dh *= tmh->de;
    tmh->mh[i] = hm + .5*dh;
    hm += dh;
  }
  return 0;
}

/* iteratively compute the density of states */
int tmh_calcdos(tmh_t *tmh, int itmax, double tol,
    const char *fndos, const char *fnlnz)
{
  int i, j, it, ie0, ie1, en = tmh->en, tpn = tmh->tpn;
  double x, dif, lnz0, db;
  double *lnn, *lnm, *bet, *lnz1;
  const double LOG0 = -1e10;
  FILE *fp;

  if (itmax <= 0) itmax = 1000;

  /* determine nonempty energy range for more efficient loops */
  for (ie0 = 0; ie0 < en; ie0++) {
    for (j = 0; j < tmh->tpn; j++)
      if (tmh->tpehis[j*en + ie0] > 0.) break;
    if (j < tmh->tpn) break;
  }
  for (ie1 = en; ie1 > ie0; ie1--) {
    for (j = 0; j < tmh->tpn; j++) 
      if (tmh->tpehis[j*en + ie1-1] > 0.) break;
    if (j < tmh->tpn) break;
  }

  /* n[j] is the total number of visits to temperature j 
   * m[i] is the total number of visits to energy i */
  xnew(lnn, tpn);
  xnew(lnm, en);
  for (i = 0; i < en; i++) lnm[i] = 0.;
  for (j = 0; j < tpn; j++) {
    lnn[j] = 0.;
    for (i = ie0; i < ie1; i++) {
      x = tmh->tpehis[j*en + i];
      lnn[j] += x;
      lnm[i] += x;
    }
  }
  for (j = 0; j < tpn; j++) 
    lnn[j] = (lnn[j] > 0.) ? log(lnn[j]) : LOG0;
  for (i = 0; i < en; i++)
    lnm[i] = (lnm[i] > 0.) ? log(lnm[i]) : LOG0;

  xnew(bet, tpn);
  for (j = 0; j < tpn; j++) 
    bet[j] = 1.0/(tmh->tp0 + (j+.5)*tmh->dtp);

  /* get mh and lnz */
  for (i = 0; i < en; i++) tmh->lng[i] = LOG0;
  tmh_calcmh(tmh);
  /* estimate initial lnz */
  for (lnz0 = 0., j = 0; j < tpn; j++) {
    x = tmh->erg0 + (j+.5)*tmh->dtp * tmh->dergdt;
    i = (int)((x - tmh->emin)/tmh->de);
    die_if (i < 0 || i >= en, "i %d, x %g\n", i, x);
    db = 1.0/(tmh->tp0 + (j+1)*tmh->dtp)
       - 1.0/(tmh->tp0 + j*tmh->dtp);
    x = tmh->mh[i]*db;
    tmh->lnz[j] = lnz0 - x*.5;
    lnz0 -= x;
  }
  for (j = tpn - 1; j >= 0; j--) 
     tmh->lnz[j] -= tmh->lnz[0];

  xnew(lnz1, tpn);

  /* repeat until convergence */
  for (it = 0; it < itmax; it++) {
    /* compute the density of states */
    for (i = ie0; i < ie1; i++) {
      for (x = LOG0, j = 0; j < tpn; j++) {
        x = lnadd(x, lnn[j] - bet[j]*tmh->mh[i] - tmh->lnz[j]);
      }
      tmh->lng[i] = (lnm[i] < LOG0+.1) ? LOG0 : (lnm[i] - x);
    }
    
    /* update partition function */
    for (j = 0; j < tpn; j++) {
      for (x = LOG0, i = ie0; i < ie1; i++)
        x = lnadd(x, tmh->lng[i] - bet[j]*tmh->mh[i]);
      lnz1[j] = x;
    }
    for (j = tpn - 1; j >= 0; j--) 
      lnz1[j] -= lnz1[0];

    /* check difference */
    for (dif = 0., j = 1; j < tpn; j++) {
      x = fabs(tmh->lnz[j] - lnz1[j]);
      if (x > dif) dif = x;
    }
    for (j = 0; j < tpn; j++) tmh->lnz[j] = lnz1[j];
    if (dif < tol) break;
  }
  
  /* write dos */
  if (fndos && (fp = fopen(fndos, "w")) != NULL) {
    for (i = ie0; i < ie1; i++) {
      x = tmh->emin + (i + .5)*tmh->de;
      fprintf(fp, "%g %g %g %g\n", x, tmh->lng[i] - tmh->lng[ie0], 
          exp(lnm[i]), tmh->mh[i]);
    }
    fclose(fp);
  }

  /* write lnz */
  if (fnlnz && (fp = fopen(fnlnz, "w")) != NULL) {
    for (j = 0; j < tpn; j++) {
      fprintf(fp, "%g %g %g\n", bet[j], tmh->lnz[j], exp(lnn[j]));
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

