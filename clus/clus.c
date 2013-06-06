#include "util.h"
#include "rng.h"
#include "mds.c"
#ifndef CLUS_C__
#define CLUS_C__
/* abstract function for cluster analysis using Monte Carlo sampling
 * given the distance matrix, we give out a cluster*/
#include "clus.h"

/* allocate space for the distance matrix */
static float **dismat_alloc(int n)
{
  float **mat, *larr;
  int npr, offset, i;

  xnew(mat, n);
  npr = (n+1)*n/2;
  xnew(larr, npr);
  for (offset = 0, i = 0; i < n; i++) {
    mat[i] = larr + offset - i;
    mat[i][i] = 0.;
    offset += n - i;
  }
  return mat;
}

/* free the distance matrix */
static void dismat_free(float **mat)
{
  free(mat[0]);
  free(mat);
}

int CLIDBLOCKSIZ = 16; /* block size for allocating clus.idx */

/* initialize a cluster of a single point i */
static void cl_init(clus_t *cl, int i, double wt)
{
  cl->cnt = 1;
  cl->cap = CLIDBLOCKSIZ;
  xnew(cl->idx, cl->cap);
  cl->idx[0] = i;
  cl->smdis = 0.;
  cl->smwt  = wt;
}

static void cl_free(clus_t *cl)
{
  if (cl->cnt > 0) {
    free(cl->idx);
    memset(cl, 0, sizeof(*cl));
  }
}

/* duplicate information in cluster b to cluster a */
static void cl_copy(clus_t *a, const clus_t *b)
{
  int k;

  a->cnt = b->cnt;
  a->cap = a->cnt;
  xrenew(a->idx, a->cap);
  for (k = 0; k < a->cnt; k++)
    a->idx[k] = b->idx[k];
  a->smdis = b->smdis;
  a->smwt  = b->smwt;
  a->centroid = b->centroid;
}

/* add a new point `i' into the cluster `cl'
 * update deldis and delwt */
static void cl_padd(clus_t *cl, int i, double deldis, double delwt)
{
  cl->cnt++;
  if (cl->cnt > cl->cap) {
    cl->cap += CLIDBLOCKSIZ;
    xrenew(cl->idx, cl->cap);
  }
  cl->idx[cl->cnt - 1] = i;
  cl->smdis += deldis;
  cl->smwt  += delwt;
}

/* remove the kth item in cluster 'cl', the energy change de  */
static void cl_premove(clus_t *cl, int k, double deldis, double delwt)
{
  die_if (k >= cl->cnt, "index error %d >= %d\n", k, cl->cnt);
  cl->cnt--;
  if (k < cl->cnt)
    cl->idx[k] = cl->idx[cl->cnt];
  cl->smdis += deldis;
  cl->smwt  += delwt;
}

#define cls_getwt(cls, i)  (cls->wt ? cls->wt[i] : 1.)

/* init. an cluster configuration from an nxn distance matrix */
clsys_t *cls_init(float **mat, double *wt, int n, double mu)
{
  clsys_t *cls;
  int ic;
  double wi;

  xnew(cls, 1);
  cls->np = n;
  cls->nc = n;
  cls->mat = mat;
  cls->wt  = wt;
  xnew(cls->c, cls->nc);
  cls->bet = 0.0;
  cls->wtot = 0.;
  for (ic = 0; ic < cls->nc; ic++) {
    wi = cls_getwt(cls, ic);
    cl_init(cls->c + ic, ic, wi);
    cls->wtot += wi;
  }
  cls->mu0 = mu;
  cls->muw = .5*mu*cls->wtot;
  cls->ene = cls->nc*cls->muw;
  return cls;
}

/* copy b to a, assume a is made by cls_init */
static void cls_copy(clsys_t *a, const clsys_t *b)
{
  int ic;
  a->np = b->np;
  a->nc = b->nc;
  a->mat = b->mat;
  a->wt = b->wt;
  a->bet = b->bet;
  a->wtot = b->wtot;
  a->mu0 = b->mu0;
  a->muw = b->muw;
  a->ene = b->ene;
  for (ic = 0; ic < a->nc; ic++) {
    cl_copy(a->c + ic, b->c + ic);
  }
}

/* find point ip in clusters,
 * return cluster id, *pk is ip's index in the cluster */
static int cls_find(clsys_t *cls, int ip, int *k)
{
  clus_t *ci;
  int ic;

  for (ic = 0; ic < cls->nc; ic++) {
    ci = cls->c + ic;
    for (*k = 0; *k < ci->cnt; (*k)++) {
      if (ip == ci->idx[*k])
        return ic;
    }
  }
  die_if (1, "cannot find point ip %d, np %d, nc %d\n", ip, cls->np, cls->nc);
  return -1;
}

/* add a new cluster with a single point */
static void cls_cadd(clsys_t *cls, int i)
{
  double wt;
  ++cls->nc;
  die_if (cls->nc > cls->np, "too many clusters %d > %d\n", cls->nc, cls->np);
  wt = cls_getwt(cls, i);
  cl_init(cls->c + cls->nc - 1, i, wt);
}

/* remove cluster ic */
static void cls_cremove(clsys_t *cls, int ic)
{
  die_if (ic >= cls->nc, "invalid cluster %d >= %d\n", ic, cls->nc);
  cls->nc--;
  if (ic < cls->nc)
    cl_copy(cls->c + ic, cls->c + cls->nc);
  cl_free(cls->c + cls->nc);
}

void cls_free(clsys_t *cls, int matwt)
{
  int ic;
  for (ic = 0; ic < cls->np; ic++)
    cl_free(cls->c + ic);
  free(cls->c);
  if (matwt) {
    if (cls->wt != NULL) free(cls->wt);
    dismat_free(cls->mat);
  }
  memset(cls, 0, sizeof(*cls));
  free(cls);
}

/* compute distance between two points i and j */
INLINE double cls_pdis(clsys_t *cls, int i, int j)
{
  die_if (i == j || i < 0 || i >= cls->np || j < 0 || j >= cls->np,
      "bad index  i: %d, j: %d, np: %g", i, j, cls->np);
  if (i < j) return cls->mat[i][j]; else return cls->mat[j][i];
}

/* check cluster configuration */
static void cls_check(clsys_t *cls, int acctype)
{
  clus_t *cl;
  int *in, ic, k, ip;

  xnew(in, cls->np);
  for (ic = 0; ic < cls->nc; ic++) {
    cl = cls->c + ic;
    die_if (cl->cnt <= 0, "invalid size %d, %d\n", ic, cl->cnt);
    for (k = 0; k < cl->cnt; k++) {
      ip = cl->idx[k];
      die_if (ip >= cls->np || in[ip], "type %d/%d: invalid ip %d/%d", acctype, cls->iter, ip, cls->np);
      in[ip] = 1;
    }
  }
  free(in);
}

/* comparison of two integers */
static int cls_cmpint_(const void *a, const void *b)
{
  return (*(const int *) a) - (*(const int *) b);
}

/* sort 1) indices in a cluster, 2) clusters by size */
static void cls_sort(clsys_t *cls)
{
  int ic, jc, im;
  clus_t *cl, cbuf;
  double wm;

  /* sort indices */
  for (ic = 0; ic < cls->nc; ic++) {
    cl = cls->c + ic;
    qsort(cl->idx, cl->cnt, sizeof(int), &cls_cmpint_);
  }
  /* sort clusters by size */
  for (ic = 0; ic < cls->nc - 1; ic++) {
    wm = cls->c[ic].smwt;
    for (im = ic, jc = ic + 1; jc  < cls->nc; jc++) {
      if (cls->c[jc].smwt > wm) {
        im = jc;
        wm = cls->c[jc].smwt;
      }
    }
    if (im != ic) {
      memcpy(&cbuf, cls->c + ic, sizeof(clus_t));
      memcpy(cls->c + ic, cls->c + im, sizeof(clus_t));
      memcpy(cls->c + im, &cbuf, sizeof(clus_t));
    }
  }
}

/* compute centroid */
static void cls_centroid(clsys_t *cls)
{
  int ic, i, ip, j, jp, imin;
  double dsm, wj, dmin;
  clus_t *cl;

  for (ic = 0; ic < cls->nc; ic++) {
    /* compute the centroid of cluster cl
     * closest to all other points */
    cl = cls->c + ic;
    if (cl->cnt <= 1) {
      cl->centroid = 0;
      continue;
    }
    for (dmin = 1e9, imin = -1, i = 0; i < cl->cnt; i++) {
      ip = cl->idx[i];
      /* compute distance from others */
      dsm = 0.;
      for (j = 0; j < cl->cnt; j++) {
        if (i == j) continue;
        jp = cl->idx[j];
        wj = cls_getwt(cls, jp);
        dsm += wj * cls_pdis(cls, ip, jp);
      }
      if (dsm < dmin) {
        imin = i;
        dmin = dsm;
      }
    }
    cl->centroid = imin;
  }
}

/* energy of a single cluster, i.e., its average distance */
static double cls_ene1(clsys_t *cls, int ic, unsigned flags)
{
  clus_t *cl = cls->c + ic;
  int k1, k2, ip, jp;
  double dis, ene, wi, wj, wtot;

  die_if(cl->cnt <= 0, "invalid cluster, n(%d) = %d\n", ic, cl->cnt);
  if (cl->cnt == 1) {
    cl->smdis = 0;
    cl->smwt = cls_getwt(cls, cl->idx[0]);
    return 0.;
  }
  for (dis = 0., wtot = 0., k1 = 0; k1 < cl->cnt; k1++) {
    ip = cl->idx[k1];
    wi = cls_getwt(cls, ip);
    wtot += wi;
    for (k2 = k1+1; k2 < cl->cnt; k2++) {
      jp = cl->idx[k2];
      wj = cls_getwt(cls, jp);
      die_if (ip == jp, "same ip jp %d\n", ip);
      dis += wi*wj*cls_pdis(cls, ip, jp);
      if (flags & CLUS_VVERBOSE)
        printf("clus %d: ip %d %g jp %d %g dis %g\n", ic, ip, wi, jp, wj, dis);
    }
  }
  if (flags & CLUS_CHECK) {
    die_if (fabs(cl->smdis - dis) > 1e-2,
      "corruption distance sum %d, %g, vs. %g\n", ic, cl->smdis, dis);
  }
  cl->smdis = dis;
  cl->smwt  = wtot;
  ene = dis/wtot;
  if (flags & CLUS_VERBOSE)
    printf("ENE: clus %d: %d points, E %g = %g/%g\n", ic, cl->cnt, ene, dis, wtot);
  return ene;
}

/* energy of all clusters */
static double cls_ene(clsys_t *cls, unsigned flags)
{
  int ic;
  double ene;

  ene = cls->nc * cls->muw;
  for (ic = 0; ic < cls->nc; ic++) {
    ene += cls_ene1(cls, ic, flags);
  }
  if (flags & CLUS_VERBOSE)
    printf("ENE: %d clusters, ene = %g\n", cls->nc, ene);
  return ene;
}

/* change mu of the cluster  */
void cls_changemu(clsys_t *cls, double mu)
{
  cls->mu0 = mu;
  cls->muw = .5 * mu * cls->wtot;
  cls->ene = cls_ene(cls, 0); /* NOTE: to improve */
}

/* merge a single-point cluster i to cluster j */
static void cls_merge(clsys_t *cls, int i, int j, double disj, double dene)
{
  clus_t *ci = cls->c + i;
  int ip;
  double wi;

/*
  printf("cls_merge, consistency check, ene %g vs %g\n", cls->ene, cls_ene(cls, 0));
*/
  die_if (ci->cnt != 1, "cluster %d is not alone %d\n", i, ci->cnt);
  ip = ci->idx[0];
  wi = cls_getwt(cls, ip);
  cl_padd(cls->c + j, ip, disj, wi); /* add ip to cluster j */
/*
  printf("cls_merge, before removal of cluster i %d, wi %g, ene %g vs. %g\n", i, wi, cls->ene, cls_ene(cls, 0));
*/
  cls_cremove(cls, i); /* remove cluster i, done after cl_padd() to avoid messing up with cluster index */
/*
  printf("cls_merge, before ene %g, dene %g, after should be %g, actually %g\n", cls->ene, dene, cls->ene + dene, cls_ene(cls, 0));
*/
  cls->ene += dene;
}

/* remove the kth point in cluster i, and add it to cluster j */
static void cls_transfer(clsys_t *cls, int i, int k, int j,
    double disi, double disj, double dene)
{
  clus_t *ci = cls->c + i;
  int ip;
  double wi;

  die_if (ci->cnt <= 1, "cluster %d is alone %d\n", i, ci->cnt);
  die_if (k < 0 || k >= ci->cnt, "point %d out of range %d\n", k, ci->cnt);
  ip = ci->idx[k];
  wi = cls_getwt(cls, ip);
  cl_premove(ci, k, -disi, -wi);
  if (j == cls->nc) {
    cls_cadd(cls, ip);
  } else {
    cl_padd(cls->c + j, ip, disj, wi);
  }
  cls->ene += dene;
}

/* remove the kth point in the cluster i,
 * and add it to the cluster j
 * combines cls_merge() and cls_transfer()
 * removes cluster i, if there is only one point there, i.e., cls_merge() */
static int cls_pmove(clsys_t *cls, int i, int ik, int j,
    double disi, double disj, double dene)
{
  const int freq = 1000;
  clus_t *ci = cls->c + i;
  int mvtype = 0, ip = ci->idx[ik];

  if (i == j) {
    return 0;
  } else if (ci->cnt == 1) {
    cls_merge(cls, i, j, disj, dene);
    mvtype = 1;
  } else {
    cls_transfer(cls, i, ik, j, disi, disj, dene);
    mvtype = 2;
  }

  cls_check(cls, mvtype);
  if (++cls->acc % freq == 0) {
    double ene1 = cls_ene(cls, CLUS_CHECK);
    die_if (fabs(ene1 - cls->ene) > 1e-1,
        "iter %d, type %d: ene diff %g, %g; de %g, clus i: %d (%d:%d), clus j: %d\n",
        cls->iter, mvtype, ene1, cls->ene, dene, i, ik, ip, j);
    cls->ene = ene1;
  }
  return mvtype;
}

/* compute the energy change of adding from point k of cluster i to cluster j */
static double cls_deadd(clsys_t *cls, int i, int k, int j, double *dis)
{
  clus_t *ci = cls->c + i, *cj = cls->c + j;
  double wi, wj;
  int ip, k1, jp;

  *dis = 0.;
  if (j == cls->nc) return cls->muw;
  ip = ci->idx[k];
  wi = cls_getwt(cls, ip);
  for (k1 = 0; k1 < cj->cnt; k1++) {
    jp = cj->idx[k1];
    wj = cls_getwt(cls, jp);
    *dis += wj*cls_pdis(cls, ip, jp);
  }
  *dis *= wi;
/*
  printf("add %d from cluster %d to %d (cnt %d, first %d): smdis %g, ddis %g, wt %g, wi %g\n", ip, i, j, cj->cnt, cj->idx[0], cj->smdis, *dis, cj->smwt, wi);
*/
  return (cj->smdis + *dis)/(cj->smwt + wi) - cj->smdis/cj->smwt;
}

/* compute the energy change of removing kth point from cluster i */
static double cls_deremove(clsys_t *cls, int i, int k, double *dis)
{
  clus_t *ci = cls->c + i;
  double wi, wj;
  int ip, k1, jp;

  die_if (i >= cls->nc, "i: %d, nc %d\n", i, cls->nc);
  if (ci->cnt == 1) {
    *dis = 0;
    return -cls->muw;
  }
  die_if (k < 0 || k >= ci->cnt, "bad index k = %d, cnt = %d\n", k, ci->cnt);
  ip = ci->idx[k];
  wi = cls_getwt(cls, ip);
  for (*dis = 0., k1 = 0; k1 < ci->cnt; k1++) {
    if (k == k1) continue;
    jp = ci->idx[k1];
    wj = cls_getwt(cls, jp);
    *dis += wj * cls_pdis(cls, ip, jp);
  }
  *dis *= wi;
  return (ci->smdis - *dis)/(ci->smwt - wi) - ci->smdis/ci->smwt;
}

/* one step, metropolis move of a cluster
 * iter is the index of iteration */
static int cls_metropolis(clsys_t *cls, unsigned flags)
{
  clus_t *ci;
  int i, j, k, acc;
  double dene, disi, disj;

  if (cls->nc < 1) return 0;
  i = (int)(cls->nc * rnd0());
  ci = cls->c + i;
  k = (int)(ci->cnt * rnd0());
  die_if(ci->cnt <= 0, "empty cluster %d\n", i);
  if (ci->cnt == 1) { /* a cluster of a single point */
    /* try to merge to another cluster */
    if ((flags & CLUS_NOGIANT) && cls->nc == 2) return 0; /* forbid giant formation */
    j = (int)((cls->nc - 1) * rnd0());
    j = (i + 1 + j) % cls->nc; /* choose any other cluster */
  } else {
    j = (int)(cls->nc * rnd0());
    j = (i + 1 + j) % (cls->nc + 1); /* m for a new cluster */
  }
  dene = cls_deremove(cls, i, k, &disi) + cls_deadd(cls, i, k, j, &disj);
/*
  printf("metro: %d, i %d, j %d, ene %g, %g; dene %g\n", cls->iter, i, j, cls_ene(cls, 0), cls->ene, dene);
*/
  acc = ((dene < 0.0) || rnd0() < exp(-cls->bet*dene));
  if ((flags & CLUS_NOGIANT) && cls->nc == 1) /* always encourage forming two clusters */
    acc = 1;
  if (acc) { /* accept the move */
    cls_pmove(cls, i, k, j, disi, disj, dene);
  }
  return 0;
}

/* choose one from a heat bath choice */
static int heatbathchoose(double *dene, int n, double bet)
{
  double *prob, demin, deb, r;
  int j;

  /* select a cluster to join */
  xnew(prob, n);
  for (demin = 1e9, j = 0; j < n; j++) {
    if (dene[j] < demin) demin = dene[j];
  }
  for (j = 0; j < n; j++) { /* probability */
    deb = bet*(dene[j] - demin);
    if (deb > 100.0) prob[j] = 0.;
    else prob[j] = exp(-deb);
  }
  for (j = 1; j < n; j++) {
    prob[j] = prob[j-1] + prob[j];
  }
  r = prob[n - 1]*rnd0();
  for (j = 0; j < n; j++) {
    if (r < prob[j]) break;
  }
  die_if (j >= n, "bad index %d > %d\n", j, n);
  free(prob);
  return j;
}

/* heat bath algorithm for ip, if ip == -1, randomly pick a particle */
static int cls_heatbath(clsys_t *cls, int ip, unsigned flags)
{
  clus_t *ci;
  int i, j, k = -1, mvtype = 0;
  double dei, disi;
  double *dene, *disj;

  if (ip >= 0 && ip < cls->np) {
    i = cls_find(cls, ip, &k);
    ci = cls->c + i;
  } else {
    i = (int)(cls->nc*rnd0());
    ci = cls->c + i;
    k = (int)(ci->cnt*rnd0());
  }
  die_if(ci->cnt <= 0 || k < 0, "empty cluster %d or bad k %d", ci->cnt, k);

  if (flags & CLUS_NOGIANT) {
    /* forbid the transition from cluster i to the other (giant) */
    if (cls->nc == 2 && ci->cnt == 1)
      return 0;
  }

  xnew(dene, cls->nc + 1); /* energy change vector */
  xnew(disj, cls->nc + 1);
  /* energy of removing k'th point from cluster i */
  dei = cls_deremove(cls, i, k, &disi);
  /* energy of adding that point to cluster j */
  for (j = 0; j <= cls->nc; j++) {
    if (j == i) continue;
    dene[j] = dei + cls_deadd(cls, i, k, j, &disj[j]);
  }
  if (ci->cnt == 1) dene[cls->nc] = 1e9; /* disable the last option for an alone cluster */
  dene[i] = 0.0;

  /* choose from the heat bath */
  if ((flags & CLUS_NOGIANT) && cls->nc == 1) {
    j = 1;
  } else {
    j = heatbathchoose(dene, cls->nc + 1, cls->bet);
  }
  if (j == i) goto EXIT; /* no move */

  /* accept the move */
  mvtype = cls_pmove(cls, i, k, j, disi, disj[j], dene[j]);
EXIT:
  free(dene);
  free(disj);
  return mvtype;
}

/* energy minimization */
static int cls_minimize(clsys_t *cls, unsigned flags)
{
  int iter, ip, changed = 1, nc0 = cls->nc, nc1, verbose = flags & (CLUS_VERBOSE|CLUS_VVERBOSE);
  double bet, ene0 = cls->ene;

  bet = cls->bet;
  cls->bet = 1e9;
  iter = cls->iter;
  for (cls->iter = 0; ; cls->iter++) {
    changed = 0;
    nc1 = cls->nc;
    for (ip = 0; ip < cls->np; ip++) {
      if (cls_heatbath(cls, ip, 0))
        changed++;
    }
    if (verbose >= 2) printf("%d: nc: %d -> %d, changed = %d\n", cls->iter, nc1, cls->nc, changed);
    if (!changed) break;
  }
  cls->bet = bet;
  if (verbose)
    printf("enemin %4d iter: (%3d %12.7f) --> (%3d %12.7f)\n",
      cls->iter, nc0, ene0, cls->nc, cls->ene);
  ip = cls->iter > 1;
  cls->iter = iter;
  return ip;
}

/* do multidimensional scaling, cls_centroid should be called
 * assuming clusters are in descending order */
static int cls_mdscal(clsys_t *cls, int nmax, int cntmin)
{
  int i, j, ip, jp, n;
  clus_t *ci, *cj;
  real *dm, *xy;

  /* construct cluster-cluster distance matrix */
  n = intmin(cls->nc, nmax);
  for (i = 0; i < n; i++) {
    if (cls->c[i].cnt <= cntmin)
      break;
  }
  if (i < n) n = i;
  if (n < 2) return 0;

  printf("doing multidimensional scaling n %d\n", n);

  xnew(dm, n*n);
  for (i = 0; i < n; i++) {
    ci = cls->c + i;
    ip = ci->idx[ci->centroid];
    for (j = i+1; j < n; j++) {
      cj = cls->c + j;
      jp = cj->idx[cj->centroid];
      dm[i*n+j] = dm[j*n+i] = (real) cls_pdis(cls, ip, jp);
    }
    dm[i*n+i] = 0.;
  }

  /* initialize coordinates */
  xnew(xy, n*2);

  /* low level multidimensional scaling */
  mds_min0(xy, dm, n, 2, 1e-16);

  /* copy coordinates */
  for (i = 0; i < n; i++) {
    ci = cls->c + i;
    ci->x = xy[i*2];
    ci->y = xy[i*2+1];
  }

  free(xy);
  free(dm);
  return 0;
}

/* prepare cls for a better representation */
static void cls_trim(clsys_t *cls)
{
  int mds_nmax = 100, cmin = 0;

  cls_sort(cls); /* sort the cluster list */
  cls_centroid(cls);
  cls_mdscal(cls, mds_nmax, cmin);
}

/* evenly spaced nc clusters */
static int cls_evensplit(clsys_t *cls, int nc)
{
  int i, ic, ppc;
  clus_t *ci;

  die_if (nc > cls->np, "too many clusters %d > %d\n", nc, cls->np);
  cls->nc = nc;

  /* compute number of points in each cluster */
  ppc = cls->np / nc;

  for (ic = 0; ic < nc; ic++) {
    ci = cls->c + ic;
    ci->cnt = (ic == nc - 1) ? (cls->np - ppc*(nc-1)) : ppc;
    ci->cap = ci->cnt;
    die_if(ci->cnt <= 0, "cnt %d, np %d, ic %d, nc %d, ppc %d\n",
        ci->cnt, cls->np, ic, nc, ppc);
    xrenew(ci->idx, ci->cap);
    for (i = 0; i < ci->cnt; i++) {
      ci->idx[i] = ic*ppc + i;
    }
  }
  for (ic = nc; ic < cls->np; ic++)
    cl_free(cls->c + ic);
  cls->ene = cls_ene(cls, 0);
  return 0;
}


/* main function for clustering: given a distance matrix return a clustering conf. */
clsys_t *cls_anneal(clsys_t *cls,
    int itmax, int method, double bet0, double bet1)
{
  for (cls->iter = 0; cls->iter < itmax; cls->iter++) {
    cls->bet = bet0 + (bet1 - bet0)*cls->iter/itmax;
    if (method == CLUS_HEATBATH) {
      cls_heatbath(cls, -1, 0);
    } else {
      cls_metropolis(cls, 0);
    }
  }
  cls_minimize(cls, CLUS_VERBOSE);
  cls_trim(cls);
  return cls;
}

clsys_t *cls_zalgo(clsys_t *cls, int itermax, int method,
    double bet0, double bet1, int nbet,
    int nstmin, int verbose)
{
  int nstage = 30, ncm = 10;
  int isz, is, ib, jb, i, it, freemode = 0;
  double lnf, lnfree, r, gemin, runmin = 1e9;
  double *barr, *lnz, *hist, *esum;
  double *muhist;
  clsys_t **clsb, *clsm, *clstmp;

  /* size dependent cluster configuration */
  if (ncm > cls->np-1) ncm = cls->np-1;
  gemin = 1e9;
  xnew(clsb, ncm+1);
  xnew(muhist, ncm+1);
  for (i = 0; i <= ncm; i++) {
    clsb[i] = cls_init(cls->mat, cls->wt, cls->np, cls->mu0);
    if (i > 0) {
      cls_evensplit(clsb[i], i);
    }
  } /* 0 is for all cluster config. w/ size > ncm */
  xnew(clsm, 1);
  clsm = cls_init(cls->mat, cls->wt, cls->np, cls->mu0);
  xnew(clstmp, 1);
  clstmp = cls_init(cls->mat, cls->wt, cls->np, cls->mu0);

  if (nbet < 2) nbet = 2;
  xnew(lnz, nbet);
  xnew(barr, nbet);
  xnew(hist, nbet);
  xnew(esum,  nbet);
  for (ib = 0; ib < nbet; ib++) {
    barr[ib] = exp(log(bet0) + ib*(log(bet1)  - log(bet0))/(nbet-1));
    lnz[ib] = 0.;
  }
  ib = 0;
  cls->bet = barr[0];
  cls->iter = 1;
  if (itermax <= 0) { /* automatically determine */
    itermax = (cls->np < 10000) ? 2*cls->np*cls->np : 10000000;
    if (verbose)
      printf("automatically determine itermax = %d for n = %d\n", itermax, cls->np);
  }
  lnfree = 0.;
  gemin = 1e9;
  for (lnf = 1., is = 0; is <= nstage; is++, lnf *= .316227766) {
    /* clear histogram */
    for (i = 0; i < nbet; i++) esum[i] = hist[i] = 0.;

    for (it = 1; ; it++) {
      int cl_flags = CLUS_NOGIANT; /* do not sample a single cluster */
      if (method == CLUS_HEATBATH) {
        cls_heatbath(cls, -1, cl_flags);
      } else {
        cls_metropolis(cls, cl_flags);
      }
      /* try a temperature transition */
      jb = (int)(ib + rnd0()*(nbet-1) + 1) % nbet;
      die_if (jb < 0 || jb > nbet, "bad jb = %d\n", jb);
      r = cls->ene*(barr[jb] - barr[ib]) + lnz[jb] - lnz[ib];
      if (r < 0. || rnd0() < exp(-r)) {
        ib = jb;
        cls->bet = barr[ib];
      }
      lnz[ib] += lnf;
      hist[ib] += 1.;
      esum[ib] += cls->ene;

      isz = (cls->nc > ncm) ? 0 : cls->nc;
      muhist[isz] += 1.;
      if (cls->ene < clsb[isz]->ene) {
        cls_copy(clsb[isz], cls);
        if (verbose >= 2)
          printf("found new emin[%d]  = %g     \r", isz, clsb[isz]->ene);
      }
      if (cls->ene < gemin) {
        gemin = cls->ene;
        if (verbose)
          printf("found new gemin = %g     \r", gemin);
      }
      if (++cls->iter >= itermax) break;

      if (it % 100 == 0) {
        /* clean up lnz */
        for (i = 0; i < nbet; i++) lnz[i] -= lnz[nbet-1];
        lnfree = 1.*nbet/cls->iter;
        if (lnf < lnfree && is > 3)
          freemode = 1;
        if (freemode) {
          /* once enter freemode never quit */
          lnf = lnfree;
        } else {  /* check histogram */
          double hmin = 1e9, hmax = 0.;
          for (i = 0; i < nbet; i++) {
            if (hist[i] < hmin) hmin = hist[i];
            if (hist[i] > hmax) hmax = hist[i];
          }
          if ((hmax-hmin)/(hmax+hmin) < 0.3) break; /* each bet met once */
        }
      }
      if (it % nstmin == 0) { /* try to minimize */
        cls_copy(clstmp, cls);
        cls_minimize(clstmp, 0);
        if (clstmp->ene < runmin) {
          if (runmin - clstmp->ene > 1e-6)
            printf("runtime emin %g, %d clusters\n", clstmp->ene, clstmp->nc);
          runmin = clstmp->ene;
          cls_copy(clsm, clstmp);
        }
      }
    }
    if (verbose) {
      printf("stage %d is complete after %d/%d, emin = %g, lnf = %g, lnfree = %g\n",
          is, it, cls->iter, gemin, lnf, lnfree);
    }
    if (verbose >= 3 || (verbose >= 2 && cls->iter >= itermax)) {
      for (i = 0; i < nbet; i++) {
        double eav;
        if (hist[i] > 0.) eav = esum[i]/hist[i]; else eav = 0.;
        printf("%8.4f %8.2f %g %g\n", barr[i], lnz[i], hist[i], eav);
      }
    }
    if (cls->iter >= itermax) {
      printf("exceeds the maximal iterations %d in stage %d\n", cls->iter, is);
      break;
    }
  }

  free(lnz);
  free(barr);
  free(hist);
  free(esum);

  /* search for the best */
  for (gemin = 1e9, i = -1, isz = 0; isz <= ncm; isz++) {
    printf("%4d, %10.0f ", isz, muhist[isz]);
    cls_minimize(clsb[isz], CLUS_VERBOSE);
    if (clsb[isz]->ene < gemin) {
      gemin = clsb[isz]->ene;
      i = isz;
    }
  }
  if (runmin < gemin) {
    gemin = runmin;
    cls_copy(cls, clsm);
  } else {
    cls_copy(cls, clsb[i]);
  }
  for (isz = 0; isz <= ncm; isz++)
    cls_free(clsb[isz], 0);
  cls_free(clsm, 0);
  free(clsb);
  free(muhist);
  cls_trim(cls);
  return cls;
}

/* write cluster results
 * `whead' is a call-back function for writing extra information */
int cls_write(clsys_t *cls, const char *fn,
    void (*whead)(FILE *, const clsys_t *cls, const void *data), const void *data,
    int version)
{
  FILE *fp;
  int ic, ni, k, ip, k1;
  clus_t *ci;
  double wtot;

  xfopen(fp, fn, "w", return -1);

  /* basic information */
  fprintf(fp, "# %d %d %g %d\n", cls->np, cls->nc, cls->mu0, version);

  /* call the callback function */
  if (whead != NULL) {
    (*whead)(fp, cls, data);
  } else {
    for (k = 0; k < cls->np; k++) {
      fprintf(fp, "# %d %g\n", k, cls_getwt(cls, k));
    }
  }

  /* write the rmsd matrix */
  for (k = 0; k < cls->np - 1; k++) {
    fprintf(fp, "# %d ", k);
    for (k1 = k+1; k1 < cls->np; k1++) {
      fprintf(fp, "%.3f ", cls->mat[k][k1]);
    }
    fprintf(fp, "\n");
  }

  for (ic = 0; ic < cls->nc; ic++) {
    ci = cls->c + ic;
    ni = ci->cnt;
    /* compute cluster weight */
    for (wtot = 0., k = 0; k < ni; k++) {
      ip = ci->idx[k];
      wtot += cls_getwt(cls, ip);
    }
    fprintf(fp, "%d %d %g %d %g %g: ",
        ic, ni, wtot, ci->idx[ci->centroid], ci->x, ci->y);
    for (k = 0; k < ni; ) {
      ip = ci->idx[k];
      for (k1 = k; k1 < ni-1; k1++) {
        if (ci->idx[k1]+1 != ci->idx[k1+1])
          break;
      }
      if (k1 == k || k == ni-1) {
        fprintf(fp, "%d ", ip);;
      } else {
        fprintf(fp, "%d-%d ", ip, ci->idx[k1]);
      }
      k = k1 + 1;
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* read a cluster file and construct a rmsd matrix */
clsys_t *cls_read(const char *fn,
    void (*rhead)(FILE *, clsys_t *, void *data), void *data)
{
  FILE *fp;
  int ic, ni, k, ip, k1, np, nc, itmp, ret, ver;
  clsys_t *cls;
  clus_t *cl;
  char buf[1024] = "", word[128], *p;
  float **dismat;
  double *wt, wtot, mu, x, y;

  xfopen(fp, fn, "r", return NULL);

  /* read in basic information */
  die_if (NULL == fgets(buf, sizeof buf, fp), "%s: no first line\n", fn);
  die_if (4 != sscanf(buf, "#%d%d%lf%d", &np, &nc, &mu, &ver),
      "%s: first line broken\n", fn);

  /* allocate space */
  xnew(wt, np);
  dismat = dismat_alloc(np);

  /* read wt */
  for (k = 0; k < np; k++) {
    die_if (NULL == fgets(buf, sizeof buf, fp), "no weight %d\n", k);
    ret = sscanf(buf, "#%d%lf", &itmp, &wt[k]);
    die_if (2 != ret, "no weight information at k = %d, ret = %d\n", k, ret);
    die_if (k != itmp, "index mismatch %d vs. %d", itmp, k);
  }

  /* read the rmsd matrix */
  for (k = 0; k < np - 1; k++) {
    ret = fscanf(fp, " #%d", &itmp);
    die_if (1 != ret, "cannot scan index k = %d, ret = %d\n", k, ret);
    die_if (k != itmp, "wrong id %d should be %d\n", itmp, k);
    for (k1 = k+1; k1 < np; k1++) {
      die_if (1 != fscanf(fp, "%f", &dismat[k][k1]),
          "cannot read matrix %d, %d\n", k, k1);
    }
  }

  /* initialize a configuration */
  cls = cls_init(dismat, wt, np, mu);

  /* read cluster configuration */
  cls->nc = nc;
  for (ic = 0; ic < nc; ic++) {
    die_if (6 != fscanf(fp, " %d%d%lf%d%lf%lf : ", &k1, &ni, &wtot, &itmp, &x, &y),
        "cannot read tag info for cluster %d", ic);
    cl = cls->c + ic;
    cl->cnt = ni;
    cl->cap = cl->cnt;
    cl->x = x;
    cl->y = y;
    xrenew(cl->idx, cl->cap);
    /* start reading indices */
    for (k = 0; k < ni; ) {
      word[0] = '\0';
      die_if (1 != fscanf(fp, "%s", word), "no atoms left, ic = %d, k = %d", ic, k);
      p = strchr(word, '-');
      if (p != NULL) {
        *p = '\0';
        ip = atoi(word);
        itmp = atoi(p+1);
        die_if (ip >= itmp, "%d >= %d\n", ip, itmp);
        die_if (k+itmp+1-ip > ni, "%d + (%d,%d) > %d", k, ip, itmp, ni);
        for (k1 = ip; k1 <= itmp; k1++)
          cl->idx[k++] = k1;
      } else {
        cl->idx[k++] = atoi(word);
      }
    }
  }
  for (ic = nc; ic < cls->np; ic++) /* free the rest */
    cl_free(cls->c + ic);
  cls->ene = cls_ene(cls, 0);
  cls_centroid(cls);

  /* start reading additional info. callback */
  rewind(fp);
  die_if(NULL == fgets(buf, sizeof buf, fp), "%s: no first line (rescan)\n", fn);
  if (rhead != NULL) (*rhead)(fp, cls, data);

  fclose(fp);
  return cls;
}

#endif /* CLUS_H__ */
