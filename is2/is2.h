#ifndef IS2_H__
#define IS2_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
  int d, l, n;
  int M, E;
  int *s; /* 0 or 1 */
} is_t;

int     is2_em(is_t *is);
int     is2_check(is_t *is);
int     is2_load(is_t *is, const char *fname);
int     is2_save(const is_t *is, const char *fname);
double  is2_exact(is_t *is, double beta, double *eav, double *cv);
int     is2_pick(const is_t *is, int *h);
int     is2_flip(is_t *is, int id, int h);
is_t   *is2_open(int l);
void    is2_close(is_t *is);

/* set transition probability */
#define is2_setproba(p, bet) { \
  double x_ = exp(-4. * bet); \
  p[2] = (UINT32) (4294967295. * x_); \
  p[4] = (UINT32) (4294967295. * x_*x_); }

/* faster macro version */
#ifdef  IS2_LB  /* L = 2^LB, N = L*L */
#define LB_     IS2_LB
#define L_      (1 << LB_)
#define N_      (L_ * L_)
#define LOWB_   (32 - 2*LB_)

#define IS2_GETH(is, id, h) { \
  unsigned ix, iy; \
  iy = id / L_, ix = id % L_; \
  h = is->s[id]*(is->s[iy*L_ + (ix+1)%L_] + is->s[iy*L_ + (ix+L_-1)%L_] \
               + is->s[(iy+1)%L_*L_ + ix] + is->s[(iy-1+L_)%L_*L_ + ix]); }
#define IS2_IRND(is, id)  id = rand32() >> LOWB_;
#define IS2_PICK(is, id, h) { IS2_IRND(is, id); IS2_GETH(is, id, h); }
#define IS2_ISEQ(is, id)  id = (id + 1) % N_;
#define IS2_PSEQ(is, id, h) { IS2_ISEQ(is, id); IS2_GETH(is, id, h); }

#define IS2_FLIP(is, id, h) { \
  is->M += (is->s[id] = -is->s[id]) * 2; \
  is->E += h * 2; }

#else

#define IS2_PICK(is, id, h)  id = is2_pick(is, &h)
#define IS2_FLIP(is, id, h)  is2_flip(is, id, h)
#endif /* IS2_LB */

#endif /* IS2_H__ */

