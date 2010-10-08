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
int     is2_load(is_t *is, const char *fname);
int     is2_save(const is_t *is, const char *fname);
double  is2_exact(is_t *is, double beta, double *eav, double *cv);
int     is2_pick(const is_t *is, int *h);
int     is2_flip(is_t *is, int id, int h);
is_t   *is2_open(int l);
void    is2_close(is_t *is);

/* faster macro version */
#ifdef  IS2_LB  /* L = 2^LB, N = L*L */
#define LB_     IS2_LB
#define L_      (1 << LB_)
#define N_      (L_ * L_)
#define LOWB_   (32 - 2*LB_)

unsigned long is2rnd_;
#define IS2_PICK(is, id, h) { \
  unsigned i, j; \
  id = (is2rnd_ = mtrand()) >> LOWB_; \
  i = id / L_, j = id % L_; \
  h = is->s[i * L_ + (j+L_-1)%L_] + is->s[i * L_ + (j+1)%L_] \
    + is->s[(i-1+L_)%L_ * L_ + j] + is->s[(i+1)%L_ * L_ + j] - 2 ; \
  if (!is->s[id]) h = -h; }

#define IS2_FLIP(is, id, h) { \
  is->M += (is->s[id] = !is->s[id]) ? 2 : -2; \
  is->E += h * 4; }

#if (LB_ <= 8)  /* up to 256 x 256 */
#define IS2_RNDRES    (1./(1 << LOWB_))
/* unused random number bits, only for small systems */
#define IS2_RNDSV()   (is2rnd_ << (2*LB_)) * (1./4294967296.)
#else
#define IS2_RNDRES    (1./4294967296.)
#define IS2_RNDSV()   rnd0()
#endif

#define IS2_INFO(is, bmax) \
  if (is->l != L_) { \
    fprintf(stderr, "is->l = %d, L = %d, LB = %d\n", is->l, L_, LB_); \
    exit(1); } \
  fprintf(stderr, "beta resolution: %g\n", IS2_RNDRES*exp(4.0*bmax)/4.0);

#else

#define IS2_PICK(is, id, h)  id = is2_pick(is, &h)
#define IS2_FLIP(is, id, h)  is2_flip(is, id, h)
#define IS2_RNDSV()          rnd0()
#define IS2_INFO(is, bmax) \
  fprintf(stderr, "beta resolution: %g\n", (1./4294967296.)*exp(4.0*bmax)/4.0);
#endif /* IS2_LB */

#endif /* IS2_H__ */

