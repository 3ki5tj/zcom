#include "util.h"
#include "rng.h"
#ifndef POTTS2_H__
#define POTTS2_H__
/* two-dimensional Potts model */

typedef struct {
  int d; /* dimension */
  int q; /* number of states of each spin */
  int l, n;
  int E;  /* potential energy */
  int *M; /* M[0..q-1] number of spins in each state */
  int *s; /* s[0..n-1], each s[i] in 0..q-1 */
  /* helper vars */
  double *accprb; /* temporary accumulated probabilities, for heat bath */
  uint32_t *uproba; /* temporary probability for MC transitions */
  double *dproba;
} potts_t;



#define PT2_SETPROBA(pt, bet) { \
  double x_ = exp(-bet), prd_; \
  prd_  = x_; pt->dproba[1] = prd_; pt->uproba[1] = (uint32_t) (4294967295. * prd_); \
  prd_ *= x_; pt->dproba[2] = prd_; pt->uproba[2] = (uint32_t) (4294967295. * prd_); \
  prd_ *= x_; pt->dproba[3] = prd_; pt->uproba[3] = (uint32_t) (4294967295. * prd_); \
  prd_ *= x_; pt->dproba[4] = prd_; pt->uproba[4] = (uint32_t) (4294967295. * prd_); \
}



/* faster macros for systems with fixed (upon compiling) size
 * to use them one must define PT2_LB and PT2_Q before including
 * PT2_PICK()/PT2_PSEQ() and PT2_FLIP() */
#ifdef  PT2_LB  /* L = 2^LB, N = L*L */
#define PT2_L   (1 << PT2_LB)
#define PT2_N   (PT2_L * PT2_L)

#define PT2_GETH(pt, id, h) { \
  unsigned ix, iy; \
  for (ix = 0; ix < PT2_Q; ix++) h[ix] = 0; \
  iy = id / PT2_L, ix = id % PT2_L; \
  h[ pt->s[iy*PT2_L + (ix+1)%PT2_L]       ]++; \
  h[ pt->s[iy*PT2_L + (ix+PT2_L-1)%PT2_L] ]++; \
  h[ pt->s[(iy+1)%PT2_L*PT2_L + ix]       ]++; \
  h[ pt->s[(iy-1+PT2_L)%PT2_L*PT2_L + ix] ]++; }
#define PT2_IRND(pt, id)  id = rand32() >> (32 - 2*PT2_LB);
/* random pick */
#define PT2_PICK(pt, id, h) { PT2_IRND(pt, id); PT2_GETH(pt, id, h); }
#define PT2_ISEQ(pt, id)  id = (id + 1) % PT2_N;
/* sequential pick */
#define PT2_PSEQ(pt, id, h) { PT2_ISEQ(pt, id); PT2_GETH(pt, id, h); }

/* change spin id from so to sn (use PT2_Q instead of pt->q) */
#define PT2_NEWFACE(pt, id, so, sn) { \
  so = pt->s[id]; sn = (so + 1 + (int)(rnd0()*(PT2_Q - 1))) % PT2_Q; }

/* change spin id from so to sn according to heat bath algorithm
 * local accprb is somehow faster */
#define PT2_HEATBATH(pt, id, so, sn, h) { \
  static double accprb[PT2_Q+1] = {0.,}; double rs_; \
  so = pt->s[id]; \
  for (sn = 0; sn < PT2_Q; sn++) accprb[sn+1] = accprb[sn] + pt->dproba[4-h[sn]]; \
  for (rs_ = accprb[PT2_Q]*rnd0(), sn = 0; sn < PT2_Q; sn++) if (accprb[sn+1] > rs_) break; \
}

#define PT2_FLIP(pt, id, so, sn, h) { \
  pt->s[id] = sn; \
  pt->M[so]--; \
  pt->M[sn]++; \
  pt->E += h[so] - h[sn];  }

#else  /* non-macro version */

#define PT2_PICK(pt, id, h)  id = pt2_pick(pt, h)
#define PT2_NEWFACE(pt, id, so, sn) { \
  so = pt->s[id]; sn = (so + 1 + (int)(rnd0()*(pt->q - 1))) % (pt->q); }
#define PT2_HEATBATH(pt, id, so, sn, h) \
  pt2_heatbath(pt, id, &so, &sn, h)
#define PT2_FLIP(pt, id, so, sn, h) {so = pt->s[id]; pt2_flip(pt, id, sn, h); }

INLINE int pt2_pick(const potts_t *pt, int h[])
{
  int i, id, ix, iy, l, lm, n, nm, *p;
  int sl, sr, sd, su;

  lm = (l = pt->l) - 1;
  nm = (n = pt->n) - l;
  id = (int)(rnd0() * n);
  iy = id / l, ix = id % l;
  p = pt->s + id;
  for (i = 0; i < pt->q; i++) h[i] = 0;
  sl = ((ix != 0 ) ? *(p-1) : *(p+lm)); h[sl]++;
  sr = ((ix != lm) ? *(p+1) : *(p-lm)); h[sr]++;
  sd = ((iy != 0 ) ? *(p-l) : *(p+nm)); h[sd]++;
  su = ((iy != lm) ? *(p+l) : *(p-nm)); h[su]++;
  return id;
}



INLINE int pt2_heatbath(potts_t *pt, int id, int *so, int *sn,
    const int h[])
{
  double rs_;
  int i, mx_ = 4;
  *so = pt->s[id];
  for (i = 0; i < pt->q; i++)
    pt->accprb[i+1] = pt->accprb[i] + pt->dproba[mx_-h[i]];
  for (rs_ = pt->accprb[pt->q]*rnd0(), i = 0; i < pt->q; i++)
    if (pt->accprb[i+1] > rs_) break;
  die_if (i >= pt->q, "no suitable selection, i = %d\n", i);
  *sn = i;
  return 0;
}



/* flip site `id' to `sn', with h different neighbors */
INLINE int pt2_flip(potts_t *pt, int id, int sn, const int h[])
{
  int so = pt->s[id];
  die_if(id >= pt->n, "id %d >= n %d\n", id, pt->n);
  die_if(sn >= pt->q || sn < 0, "invalid sn %d (q = %d)\n", sn, pt->q);
  pt->s[id] = sn;
  pt->M[so]--;
  pt->M[sn]++;
  return pt->E += h[so] - h[sn];
}

#endif /* PT2_LB */



/* compute the total energy and magnetization */
INLINE int pt2_em(potts_t *pt)
{
  int i, j, l, s, s1, s2, *p;

  pt->E = 0;
  p = pt->s;
  l = pt->l;
  for (i = 0; i < pt->q; i++) pt->M[i] = 0;
  for (i = 0; i < l; i++)
    for (j = 0; j < l; j++) {
      s = p[i*l + j];
      s1 = p[((i+1)%l)*l + j];
      if (s1 == s) pt->E--;
      s2 = p[i*l + (j+1)%l];
      if (s2 == s) pt->E--;
      pt->M[s]++;
    }
  return pt->E;
}



/* check */
INLINE int pt2_check(potts_t *pt)
{
  int i, e, *mg, q = pt->q;

  for (i = 0; i < pt->n; i++) /* check spin value */
    die_if (pt->s[i] < 0 || pt->s[i] >= q, "s[%d] %d, q %d\n", i, pt->s[i], q);
  e = pt->E;
  xnew(mg, q);
  for (i = 0; i < q; i++) mg[i] = pt->M[i];
  die_if (e != pt2_em(pt), "error: E = %d, should be %d\n", e, pt->E);
  for (i = 0; i < q; i++)
    die_if (mg[i] != pt->M[i], "error: M(%d) = %d, should be %d", i, mg[i], pt->M[i]);
  free(mg);
  return 0;
}



/* pick a random site (return its id)
 * compute h[j], the numbers of neighboring spins with value j */
/* load spin configuration */
INLINE int pt2_load(potts_t *pt, const char *fname)
{
  FILE *fp;
  int i, lx, ly, n, c;
  char s[80];

  xfopen(fp, fname, "r", return -1);
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "missing first line %s\n", fname);
    return -1;
  }
  if (4 != sscanf(s, "%d%d%d%d", &i, &lx, &ly, &n)
      || i != 2 || lx != ly || lx != pt->l || n != pt->n) {
    fprintf(stderr, "bad setting: %dD, %dx%d = %d\n", i, lx, ly, n);
    return -1;
  }
  for (i = 0; i < n; i++) {
    while ((c = fgetc(fp)) != EOF && c == '\n') ;
    if (c == EOF) break;
    c -= '0';
    if (c < 0 || c >= pt->q) {
      fprintf(stderr, "BAD %s s[%d] = %d, q = %d\n", fname, i, c, pt->q);
      break;
    }
    pt->s[i] = c;
  }
  if (i < n) {
    fprintf(stderr, "%s: data stopped at i = %d, clear\n", fname, i);
    for (i = 0; i < n; i++) pt->s[i] = 0;
  }
  fclose(fp);
  pt2_em(pt); /* re-compute energy/magnetization */
  return 0;
}



/* save spin configuration */
INLINE int pt2_save(const potts_t *pt, const char *fname)
{
  FILE *fp;
  int i, j, l, *p;

  xfopen(fp, fname, "w", return -1);
  l = pt->l;
  fprintf(fp, "%d %d %d %d\n", pt->d, l, l, pt->n);
  for (p = pt->s, i = 0; i < l; i++) {
    for (j = 0; j < l; j++, p++)
      fprintf(fp, "%c", *p+'0');
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



/* initialize an lxl q-state Potts model */
INLINE potts_t *pt2_open(int l, int q)
{
  int i, n;
  potts_t *pt;

  xnew(pt, 1);
  pt->d = 2;
  pt->q = q;
  pt->l = l;
  pt->n = n = l*l;
  xnew(pt->s, n);
  xnew(pt->M, q);
  for (i = 0; i < n; i++)
    pt->s[i] = 0;
  for (pt->M[0] = n, i = 1; i < q; i++)
    pt->M[i] = 0;
  pt->E = -pt->d * n;
  xnew(pt->accprb, q+1);
  pt->accprb[0] = 0.;
  /* dynamic array of uproba/dproba seems to be faster */
  xnew(pt->uproba, 2*pt->d+1);
  pt->uproba[0] = 0xffffffffu;
  xnew(pt->dproba, 2*pt->d+1);
  pt->dproba[0] = 1.;
  return pt;
}



INLINE void pt2_close(potts_t *pt)
{
  if (pt != NULL) {
    free(pt->s);
    free(pt->M);
    free(pt->accprb);
    free(pt->uproba);
    free(pt->dproba);
    free(pt);
  }
}

#endif

