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
#define is2_setproba(p, bet) { p[2]=exp(-4.*bet); p[4]=p[2]*p[2]; }

/* faster macro version */
#ifdef  IS2_LB  /* L = 2^LB, N = L*L */
#define LB_     IS2_LB
#define L_      (1 << LB_)
#define N_      (L_ * L_)
#define LOWB_   (32 - 2*LB_)

unsigned long is2rnd_;
#define IS2_PICK(is, id, h) { \
  unsigned ix, iy; \
  id = (is2rnd_ = mtrand()) >> LOWB_; \
  iy = id / L_, ix = id % L_; \
  h = is->s[id]*(is->s[iy*L_ + (ix+1)%L_] + is->s[iy*L_ + (ix+L_-1)%L_] \
               + is->s[(iy+1)%L_*L_ + ix] + is->s[(iy-1+L_)%L_*L_ + ix]); }

#define IS2_FLIP(is, id, h) { \
  is->M += (is->s[id] = -is->s[id]) * 2; \
  is->E += h * 2; }

#if (LB_ <= 6)  /* up to 64 x 64 */
#define IS2_RNDRES    (1./(1 << LOWB_))
/* unused random number bits, only for small systems */
#define IS2_RNDSV()   (is2rnd_ << (2*LB_)) * (1./4294967296.)
#else
#define IS2_RNDRES    (1./4294967296.)
#define IS2_RNDSV()   rnd0()
#endif

#define IS2_INFO(is, beta) { double tmp = exp(4.0*beta)*.25; \
  if (is->l != L_) { \
    fprintf(stderr, "error: is->l = %d, L = %d, LB = %d\n", is->l, L_, LB_); \
    exit(1); } \
  fprintf(stderr, "beta resolution: %g, %g\n", \
      IS2_RNDRES*tmp, (1./4294967296.)*tmp); }

#else

#define IS2_PICK(is, id, h)  id = is2_pick(is, &h)
#define IS2_FLIP(is, id, h)  is2_flip(is, id, h)
#define IS2_RNDSV()          rnd0()
#define IS2_INFO(is, beta) \
  fprintf(stderr, "beta resolution: %g\n", (1./4294967296.)*exp(4.0*beta)/4.0);
#endif /* IS2_LB */

#endif /* IS2_H__ */

