/* integral identity for potential energy distribution */
#define ZCOM_PICK
#define ZCOM_DISTR
#define ZCOM_CFG
#define ZCOM_RC /* used by distr2d */
#include "zcom.h"
#include "../distr2d.h"

const char *fncfg = "bb.cfg";
const char *fnraw = "bbcompat.txt";
const char *fnds = "bbii.txt";

int n = 360;
double xmin = -M_PI, xmax = M_PI, xdel = 2*M_PI/360;

int halfwin = 2; /* half window size */
int iitype = 0; /* 0: Adib-Jarsynski; 1: modulated */
int mfhalfwin = 0; /* half window size for mean force */

static void doii(distr2d_t *d)
{
  double mfsig;

  mfsig = sqrt( distr2d_mfvar(d) );
  printf("mfsig = %g, window size %g deg\n", mfsig, 180/M_PI/mfsig);
  distr2d_mfwin(d, 2);
  distr2d_calclnrho(d);
  distr2d_winfixed(d, d->win, 2, 1.0);
  distr2d_ii0(d, d->win);
}

int main(void)
{
  distr2d_t *d;

  d = distr2d_open(xmin, xmax, xdel, xmin, xmax, xdel);
  die_if(0 != distr2d_load(d, fnraw), "failed to load data from %s\n", fnraw);
  doii(d);
  printf("tot %g, norm %g\n", d->tot, d->tot * d->dx * d->dy);
  distr2d_save(d, fnds);
  distr2d_close(d);
  return 0;
}

