#include "util.c"
#include "hist.c"
#include "rng.c"

#ifndef TMH_C__
#define TMH_C__

/* tempering with modified Hamiltonian */
#include "tmh.h"

/* 0: cold; 1: hot */
tmh_t *tmh_open(double tp0, double tp1, double erg0, double erg1,
    double emin, double emax, double de, double ensexp,
    int dhdeorder)
{
  tmh_t *tmh;
  int i;

  xnew(tmh, 1);
  die_if(tp0 >= tp1, "Error: T0 %g >= T1 %g\n", tp0, tp1);
  tmh->tp0 = tp0;
  tmh->tp1 = tp1;

  tmh->de = de;
  tmh->emin = dblround(emin, de);
  tmh->emax = dblround(emax, de);
  tmh->en = (int)((tmh->emax - tmh->emin)/de + .5); /* number of energy bins */

  /* the updating energy range */
  die_if(erg0 >= erg1, "Error: erg0 %g >= erg1 %g\n", erg0, erg1);
  tmh->erg0 = dblround(erg0, de);
  tmh->erg1 = dblround(erg1, de);
 
  tmh->dhdeorder = dhdeorder; 
  tmh->dhdemin = 0.1;
  tmh->dhdemax = 10.0;
  xnew(tmh->dhde, tmh->en + 1);
  for (i = 0; i <= tmh->en; i++)
    tmh->dhde[i] = 1.;
  tmh->dergdt = (tmh->erg1 - tmh->erg0)/(tmh->tp1 - tmh->tp0);
  tmh->ie0 = (int)((tmh->erg0 - tmh->emin)/de + .5);
  tmh->ie1 = (int)((tmh->erg1 - tmh->emin)/de + .5) - 1;

  tmh->tpn = tmh->ie1 - tmh->ie0 + 1;
  tmh->dtp = (tmh->tp1 - tmh->tp0)/tmh->tpn * (1. + 1e-12);
  xnew(tmh->tpehis, tmh->tpn*tmh->en);

  if (tmh->dhdeorder > 0) tmh->ie1++; /* interpolation */  

  tmh->ensexp = ensexp;
  tmh_settp(tmh, (tmh->tp0+tmh->tp1)*.5);

  return tmh;
}

void tmh_close(tmh_t *tmh)
{
  if (tmh != NULL) {
    free(tmh->dhde);
    free(tmh->tpehis);
    free(tmh);
  }
}

static double tmh_hdif0(tmh_t *tmh, double eh, double el)
{
  int ie, iel, ieh;
  double dh;

  if (eh < tmh->erg0) { /* first energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ie0];
  } else if (el > tmh->erg1) { /* last energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ie1];
  } else {
    dh = 0.;
    /* overhead index */
    if (el < tmh->erg0) {
      dh += (tmh->erg0 - el)*tmh->dhde[tmh->ie0];
      el = tmh->erg0 + 1e-8;
    }
    if (eh >= tmh->erg1) {
      dh += (eh - tmh->erg1)*tmh->dhde[tmh->ie1];
      eh = tmh->erg1 - 1e-8;
    }
    /* energy index */
    iel = (int)((el - tmh->emin)/tmh->de);
    ieh = (int)((eh - tmh->emin)/tmh->de);
    if (iel == ieh) {
      dh += (eh - el)*tmh->dhde[iel];
    } else if (iel < ieh) {
      /* dh at the two terminal energy bin */
      dh += (tmh->emin + (iel+1)*tmh->de - el) * tmh->dhde[iel]
          + (eh - (tmh->emin + tmh->de*ieh)) * tmh->dhde[ieh];
      /* integrate dH/dE */
      for (ie = iel+1; ie < ieh; ie++)
        dh += tmh->dhde[ie]*tmh->de;
    }
  }
  return dh; 
}

static double tmh_hdif1(tmh_t *tmh, double eh, double el)
{
  int ie, iel, ieh;
  double dh, de, k, el0, eh0;

  if (eh < tmh->erg0) { /* first energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ie0];
  } else if (el > tmh->erg1) { /* last energy bin */
    dh = (eh - el)*tmh->dhde[tmh->ie1];
  } else {
    dh = 0.;
    if (el < tmh->erg0) {
      dh += (tmh->erg0 - el)*tmh->dhde[tmh->ie0];
      el = tmh->erg0 + 1e-8;
    }
    if (eh >= tmh->erg1) {
      dh += (eh - tmh->erg1)*tmh->dhde[tmh->ie1];
      eh = tmh->erg1 - 1e-8;
    }
    /* energy index */
    iel = (int)((el - tmh->emin)/tmh->de);
    ieh = (int)((eh - tmh->emin)/tmh->de);
    if (iel == ieh) {
      k = (tmh->dhde[iel+1] - tmh->dhde[iel])/tmh->de;
      el0 = tmh->emin + iel*tmh->de;
      dh += (eh - el)*(tmh->dhde[iel] + k * (.5f*(el+eh) - el0));
    } else if (iel < ieh) {
      /* dh at the two terminal energy bin */
      el0 = tmh->emin + (iel+1)*tmh->de; 
      de = el0 - el;
      k = (tmh->dhde[iel+1] - tmh->dhde[iel])/tmh->de;
      dh += de * (tmh->dhde[iel+1] - k*.5*de);
      eh0 = tmh->emin + tmh->de*ieh;
      de = eh - eh0;
      k = (tmh->dhde[ieh+1] - tmh->dhde[ieh])/tmh->de;
      dh += de * (tmh->dhde[ieh] + k*.5*de);
      /* integrate dH/dE */
      for (ie = iel+1; ie < ieh; ie++)
        dh += .5*(tmh->dhde[ie] + tmh->dhde[ie+1])*tmh->de;
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
  int ie, j; 
  FILE *fp;
  double ehis;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
  for (ie = 0; ie < tmh->ie0; ie++)
    tmh->dhde[ie] = tmh->dhde[tmh->ie0];
  for (ie = tmh->ie1+1; ie <= tmh->en; ie++) 
    tmh->dhde[ie] = tmh->dhde[tmh->ie1];
  fprintf(fp, "# %d %d %d %g %g %g %g %g\n", tmh->en, tmh->ie0, 
      tmh->ie1 + ((tmh->dhdeorder == 0) ? 1 : 0), 
      tmh->emin, tmh->de, amp, t, tmh->tp);
  for (ie = 0; ie < tmh->en; ie++) {
    /* compute overall energy histogram */
    for (ehis = 0., j = 0; j < tmh->tpn; j++) 
      ehis += tmh->tpehis[j*tmh->en + ie];
    fprintf(fp, "%g %g %g\n", 
        tmh->emin + ie*tmh->de, tmh->dhde[ie], ehis);
  }
  fclose(fp);
  return 0;
}

/* read dhde and overall energy distribution */
int tmh_loaddhde(tmh_t *tmh, const char *fn, double *amp, double *t)
{
  int ie, en, ie0, ie1;
  FILE *fp;
  char s[1024];
  double emin, de, erg, erg0, dhde;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "cannot read the first line %s\n", fn);
    goto ERR;
  }
  if (8 != sscanf(s+1, "%d%d%d%lf%lf%lf%lf%lf", 
        &en, &ie0, &ie1, &emin, &de, amp, t, &tmh->tp)) {
    fprintf(stderr, "corrupted info line %s", s);
    goto ERR;
  }
  if (tmh->dhdeorder == 0) ie1--;
  if (en != tmh->en || ie0 != tmh->ie0 || ie1 != tmh->ie1
   || fabs(emin - tmh->emin) > 1e-5 || fabs(de - tmh->de) > 1e-5) {
    fprintf(stderr, "bad energy en %d %d, ie0 %d %d, ie1 %d %d, emin %g %g, de %g %g\n",
        en, tmh->en, ie0, tmh->ie0, ie1, tmh->ie1, emin, tmh->emin, de, tmh->de);
    goto ERR;
  }

  for (ie = 0; ie < tmh->en; ie++) {
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "cannot read line %d\n", ie);
      goto ERR;
    }
    if (2 != sscanf(s, "%lf%lf", &erg, &dhde)) {
      fprintf(stderr, "cannot read energy and dhde at %d\n, s = %s", ie, s);
      goto ERR;
    }
    erg0 = tmh->emin + ie*tmh->de;
    if (fabs(erg0 - erg) > tmh->de*.1) {
      fprintf(stderr, "energy %g, should be %g\n", erg, erg0);
      goto ERR;
    }
    tmh->dhde[ie] = dhde;
  }
  tmh->dhde[tmh->en] = tmh->dhde[tmh->en-1];
  fclose(fp);
  return 0;
ERR:
  fclose(fp);
  return -1;
}

/* get energy range from dhde file */
int tmh_loaderange(const char *fn, double *erg0, double *erg1, 
    double *emin, double *emax, double *de)
{
  int en, ie0, ie1;
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
  if (5 != sscanf(s+1, "%d%d%d%lf%lf", &en, &ie0, &ie1, emin, de)) {
    fprintf(stderr, "corrupted info line %s", s);
    goto ERR;
  }
  *erg0 = *emin + ie0*(*de);
  *erg1 = *emin + ie1*(*de);
  *emax = *emin + en*(*de);
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
  for (i = 0; i < tmh->tpn; i++) {
    eh = tmh->tpehis + i*tmh->en;
    for (cnt = esm = 0., j = 0; j < tmh->en; j++) {
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

#endif

