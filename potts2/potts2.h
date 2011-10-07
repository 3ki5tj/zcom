#define INLINE __inline static
#ifndef POTTS2_H__
#define POTTS2_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

int pt2_em(potts_t *pt);
int pt2_check(potts_t *pt);
int pt2_load(potts_t *pt, const char *fname);
int pt2_save(const potts_t *pt, const char *fname);
potts_t *pt2_open(int l, int q);
void pt2_close(potts_t *pt);

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


#endif

