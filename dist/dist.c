#ifndef DIST_C__
#define DIST_C__

#include "dist.h"

/* to normalize and write histograms to file
 pointer 'h' (to be converted from 1D or 2D array) contains 'rows' histograms,
 each contains 'cols' entries, from 'base' to 'base+inc*cols' */
int wdistex(double *h, int rows, int cols, double base, double inc,
    int flag, const char *fname)
{
  const char *filename;
  FILE *fp;
  int i, j, imax, imin;
  double sum, *p, delta;

  filename = (fname != NULL) ? fname : "HIST";

  if ((fp = fopen(filename, "w")) == NULL) {
    printf("cannot write history file [%s].\n", filename);
    return 1;
  }
  delta = (flag & WD_ADDAHALF) ? 0.5 : 0;

  for (j = 0; j < rows; j++) {
    p = h+j*cols;

    if (flag & WD_KEEPEDGE) {
      imin = 0;
      imax = cols;
    } else {
      for (i = cols-1; i >= 0; i--)
        if (p[i] > 0)
          break;
      imax = i+1;
      if (imax == 0)
        continue;

      for (i = 0; i < imax; i++)
        if (p[i] > 0)
          break;
      imin = i;
    }

    for (sum = 0, i = imin; i < imax; i++)
      sum += p[i];
    sum *= inc;
    if (fabs(sum) < 1e-6)
      sum = 1;

    for (i = imin; i < imax; i++) {
      if ((flag & WD_NOZEROES) && p[i] < 1e-6)
        continue;
      fprintf(fp,"%g %20.14E %d\n", base+(i+delta)*inc, p[i]/sum, j);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  return 0;
}

#endif /* ZCOM_DIST__ */

