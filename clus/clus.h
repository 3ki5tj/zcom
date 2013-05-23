#ifndef CLUS_H__
#define CLUS_H__
#include <stdio.h>
#include <stdlib.h>

typedef struct { /* structure of a single cluster */
  int *idx;
  int cnt; /* number of points in this cluster, size of idx */
  int cap; /* capacity for holding */
  int centroid; /* centroid */
  double smdis; /* sum of pair distances, wi*wj*dij */
  double smwt;  /* sum of weights, wi */
  double x, y;  /* multidimensional scaling */
} clus_t;

typedef struct { /* clsys: array of all clusters */
  int np;       /* # of points */
  int nc;       /* # of clusters */
  float **mat;  /* distance matrix i < j */
  double *wt;   /* weight */
  double wtot;  /* total weight */
  clus_t *c;    /* cluster array */
  double ene;   /* energy of all clusters */
  /* auxiliary variables */
  double mu0;   /* input mu */
  double muw;   /* penalty of adding a cluster the actual one
                   the weighted version, 0.5 * mu0  */
  double bet;   /* inverse temperature */
  int acc; /* accepted metropolis moves */
  int iter;  /* iteration index */
} clsys_t;


#define CLUS_METROPOLIS 0x00
#define CLUS_VERBOSE    0x01
#define CLUS_VVERBOSE   0x02
#define CLUS_HEATBATH   0x10
#define CLUS_CHECK      0x20
#define CLUS_NOGIANT    0x100 /* ignore a single-cluster configuration during sampling */

clsys_t *cls_init(float **mat, double *wt, int n, double mu);
void cls_free(clsys_t *cls, int matwt);
void cls_changemu(clsys_t *cls, double mu);
clsys_t *cls_read(const char *fn, void (*rhead)(FILE *, clsys_t *, void *data), void *data);
int cls_write(clsys_t *cls, const char *fn,
    void (*whead)(FILE *, const clsys_t *cls, const void *data), const void *data, int ver);
clsys_t *cls_anneal(clsys_t *cls, int itmax, int method, double bet0, double bet1);
clsys_t *cls_zalgo(clsys_t *cls, int itmax, int method,
    double bet0, double bet1, int nbet, int nstmin, int verbose);

#endif

