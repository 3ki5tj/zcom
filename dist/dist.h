#ifndef DIST_H__
#define DIST_H__

#include <stdio.h>
#include <math.h>

#define WD_ADDAHALF   0x0001
#define WD_KEEPEDGE   0x0002
#define WD_NOZEROES   0x0004
#define wdist(h,rows,cols,base,inc,fname) \
  wdistex((double *)h,rows,cols,base,inc,WD_ADDAHALF,fname)

int wdistex(double *h, int rows, int cols, double base, double inc, 
    int flags, const char *fname);

#endif

