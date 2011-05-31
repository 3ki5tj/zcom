#ifndef HIST_H__
#define HIST_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define HIST_VERBOSE    0x0001

#define HIST_ADDAHALF   0x0010
#define HIST_NOZEROES   0x0020
#define HIST_KEEPLEFT   0x0040
#define HIST_KEEPRIGHT  0x0080
#define HIST_KEEPEDGE   (HIST_KEEPLEFT|HIST_KEEPRIGHT)
#define HIST_KEEPHIST  0x0100

#define HIST_ADDITION   0x1000

/* old names */
#define wdist(h,m,n,x0,dx,fn) wdistex(h,m,n,x0,dx,HIST_ADDAHALF|HIST_VERBOSE,fn)
#define wdistex histsave

#define histsave(h,rows,xcnt,xmin,dx,flags,fname) \
  histsavex((const double *)h,rows,xcnt,xmin,dx,flags,NULL,NULL,NULL,fname)

int histsavex(const double *h, int rows, int xcnt, double xmin, double dx, 
    unsigned flags, int (*fwheader)(FILE *, void *),
    double (*fnorm)(int, int, double, double, void *), 
    void *pdata, const char *fname);

#define histload(h,rows,xcnt,xmin,dx,flags,fname) \
  histloadx((double *)h,rows,xcnt,xmin,dx,flags,NULL,NULL,NULL,fname)

int histloadx(double *hist, int rows, int xcnt, double xmin, double dx,
    unsigned flags, int (*frheader)(const char *, void *),
    double (*fnorm)(int, int, double, double, void *), 
    void *pdata, const char *fn);

int histadd(const double *x, double w, double *h, int rows, 
    int xcnt, double xmin, double dx, unsigned flags);

/* object oriented wrapper functions */
typedef struct {
  double *arr;
  int rows;
  int xcnt;
  double xmin;
  double dx;
  int (*fwheader)(FILE *, void *);
  int (*frheader)(const char *, void *);
  double (*fnorm)(int, int, double, double, void *);
} hist_t;

#define hs_init(m,x0,x1,dx) hs_initx(m,x0,x1,dx,NULL,NULL,NULL)
#define hs_save(hs,fn,flags) hs_savex(hs,fn,NULL,flags)
#define hs_load(hs,fn,flags) hs_loadx(hs,fn,NULL,flags)

hist_t *hs_initx(int rows, double xmin, double xmax, double dx,
    int (*fwh)(FILE *, void *), int (*frh)(const char*, void *),
    double (*fnorm)(int, int, double, double, void *));
void hs_free(hist_t *hs);
int hs_savex(const hist_t *hs, const char *fname, void *pdata, unsigned flags);
int hs_loadx(hist_t *hs, const char *fname, void *pdata, unsigned flags);
int hs_add(hist_t *hs, const double *x, double w, unsigned flags);
int hs_add1(hist_t *hs, int j, double x, double w, unsigned flags);

#endif

