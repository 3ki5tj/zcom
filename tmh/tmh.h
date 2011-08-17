#ifndef TMH_H__
#define TMH_H__

typedef struct {
  double tp, ec; /* current temperature, and the expected energy there */
  int itp, iec; /* indices of tp and ec */
  double tp0, tp1, dtp; /* temperature range */
  int tpn; /* number of temperature */
  double erg0, erg1; /* energy range (erg0, erg1) */
  int ie0, ie1; /* correspond to erg0 and erg1 */
  double emin, emax, de; /* energy range */
  int en; /* number of energy bins */
  double dtderg; /* (tp1 - tp0)/(erg1 - erg0) */
  double dedt; /* (emax - emin)/(tp1 - tp0)/de */
  int dhdeorder;
  double *dhde; /* dH / dE - 1 */
  double *ehis; /* energy histogram */
  double *tphis, *tpesm; /* temperature histogram */
  double ensexp; /* w(T) = 1/T^ensexp */
} tmh_t;

tmh_t *tmh_open(double tp0, double tp1, double erg0, double erg1,
    double emin, double emax, double de, double ensexp, 
    int dhdeorder);
void tmh_close(tmh_t *tmh);
double tmh_hdif(tmh_t *tmh, double eb, double ea);
int tmh_tlgvmove(tmh_t *tmh, double enow, double lgvdt);
int tmh_writeerg(tmh_t *tmh, const char *fname);
int tmh_writetp(tmh_t *tmh, const char *fname);

/* set the current temperature */
ZCINLINE void tmh_settp(tmh_t *tmh, double tp)
{
  tmh->tp = tp;
  tmh->itp = (int)((tmh->tp - tmh->tp0)/tmh->dtp);
  tmh->ec = tmh->erg0 + (tmh->tp - tmh->tp0)/tmh->dtderg;
  tmh->iec = (int)((tmh->ec - tmh->emin)/tmh->de);
}

/* retrieve local dhde */
ZCINLINE double tmh_getdhde(tmh_t *tmh)
{
  if (tmh->dhdeorder == 0) {
    return tmh->dhde[tmh->iec];
  } else {
    double lam = (tmh->ec - (tmh->emin + tmh->iec*tmh->de))/tmh->de;
    return tmh->dhde[tmh->iec]*(1-lam) + tmh->dhde[tmh->iec+1]*lam;
  }
}

/* update dhde curve */
ZCINLINE void tmh_dhdeupdate(tmh_t *tmh, double erg, double amp)
{
  double del = amp * (erg - tmh->ec);
  if (tmh->dhdeorder == 0) {
    tmh->dhde[tmh->iec] += del;
  } else {
    tmh->dhde[tmh->iec] += .5*del;
    tmh->dhde[tmh->iec+1] += .5*del;
  }
}

ZCINLINE void tmh_eadd(tmh_t *tmh, double erg)
{
  int ie = (int)((erg - tmh->emin)/tmh->de);
  tmh->ehis[ie] += 1.;
}

ZCINLINE void tmh_tadd(tmh_t *tmh, double erg)
{
  tmh->tphis[tmh->itp] += 1.;
  tmh->tpesm[tmh->itp] += erg;
}

#endif /* TMH_H__ */

