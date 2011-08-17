#include "util.c"

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
  tmh->tpn = 100;
  tmh->dtp = (tmh->tp1 - tmh->tp0)/tmh->tpn * 1.000001;

  tmh->de = de;
  tmh->emin = dblround(emin, de);
  tmh->emax = dblround(emax, de);
  tmh->en = (int)((tmh->emax - tmh->emin)/de + .5); /* number of energy bins */
  /* set up the energy-to-beta ratio */
  tmh->dedt = (tmh->emax - tmh->emin)/(tmh->tp1 - tmh->tp0)/tmh->de;

  /* the updating energy range */
  die_if(erg0 >= erg1, "Error: erg0 %g >= erg1 %g\n", erg0, erg1);
  tmh->erg0 = dblround(erg0, de);
  tmh->erg1 = dblround(erg1, de);
 
  tmh->dhdeorder = dhdeorder; 
  xnew(tmh->dhde, tmh->en + 1);
  for (i = 0; i <= tmh->en; i++)
    tmh->dhde[i] = 1.;
  xnew(tmh->ehis, tmh->en);
  xnew(tmh->tphis, tmh->tpn);
  xnew(tmh->tpesm, tmh->tpn);
  tmh->dtderg = (tmh->tp1 - tmh->tp0)/(tmh->erg1 - tmh->erg0);
  tmh->ie0 = (int)((tmh->erg0 - tmh->emin)/de + .5);
  tmh->ie1 = (int)((tmh->erg1 - tmh->emin)/de + .5) - 1;
  if (tmh->dhdeorder > 0) tmh->ie1++; /* interpolation */  

  tmh->ensexp = ensexp;
  return tmh;
}

void tmh_close(tmh_t *tmh)
{
  if (tmh != NULL) {
    free(tmh->dhde);
    free(tmh->ehis);
    free(tmh->tphis);
    free(tmh->tpesm);
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

  die_if (e0 < tmh->emin || e1 > tmh->emax, 
      "energy (%g, %g) out of range (%g, %g) \n", 
      e0, e1, tmh->emin, tmh->emax);

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

int tmh_writeerg(tmh_t *tmh, const char *fname)
{
  int ie; 
  FILE *fp;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fname);
    return -1;
  }
  for (ie = 0; ie < tmh->ie0; ie++)
    tmh->dhde[ie] = tmh->dhde[tmh->ie0];
  for (ie = tmh->ie1+1; ie < tmh->en; ie++) 
    tmh->dhde[ie] = tmh->dhde[tmh->ie1];
  for (ie = 0; ie < tmh->en; ie++) {
    fprintf(fp, "%g %g %g\n", 
        tmh->emin + ie*tmh->de, tmh->dhde[ie], tmh->ehis[ie]);
  }
  fclose(fp);
  return 0;
}

int tmh_writetp(tmh_t *tmh, const char *fname)
{
  int i;
  double eav;
  FILE *fp;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fname);
    return -1;
  }
  for (i = 0; i < tmh->tpn; i++) {
    eav = (tmh->tphis[i] > 0) ? tmh->tpesm[i]/tmh->tphis[i] : 0.;
    fprintf(fp, "%g %g %g\n", 
        tmh->tp0 + i*tmh->dtp, tmh->tphis[i], eav);
  }
  fclose(fp);
  return 0;
}

#endif

