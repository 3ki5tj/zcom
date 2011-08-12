#include "rng.h"
#ifndef ISING2_H__
#define ISING2_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
  int d, l, n;
  int M, E;
  int *s; /* 0 or 1 */
  /* helper vars */
  uint32_t *uproba; /* temporary probability for MC transitions */ 
} ising_t;

int     is2_em(ising_t *is);
int     is2_check(ising_t *is);
int     is2_load(ising_t *is, const char *fname);
int     is2_save(const ising_t *is, const char *fname);
double  is2_exact(ising_t *is, double beta, double *eav, double *cv);
int     is2_pick(const ising_t *is, int *h);
int     is2_flip(ising_t *is, int id, int h);
ising_t*is2_open(int l);
void    is2_close(ising_t *is);

/* set transition probability */
#define IS2_SETPROBA(is, bet) { \
  double x_ = exp(-4. * bet); \
  is->uproba[2] = (uint32_t) ((double)(0xffffffff) * x_); \
  is->uproba[4] = (uint32_t) ((double)(0xffffffff) * x_*x_); }

/* faster macros for systems with fixed (upon compiling) size
 * to use them one must define IS2_LB before including 
 * IS2_PICK()/IS2_PSEQ() and IS2_FLIP() */
#ifdef  IS2_LB  /* L = 2^LB, N = L*L */
#define IS2_L   (1 << IS2_LB)
#define IS2_N   (IS2_L * IS2_L)

#define IS2_GETH(is, id, h) { \
  unsigned ix, iy; \
  iy = id / IS2_L, ix = id % IS2_L; \
  h = is->s[id]*(is->s[iy*IS2_L + (ix+1)%IS2_L] \
               + is->s[iy*IS2_L + (ix+IS2_L-1)%IS2_L] \
               + is->s[(iy+1)%IS2_L*IS2_L + ix] \
               + is->s[(iy-1+IS2_L)%IS2_L*IS2_L + ix]); }
#define IS2_IRND(is, id)  id = rand32() >> (32 - 2*IS2_LB);
/* random picking */
#define IS2_PICK(is, id, h) { IS2_IRND(is, id); IS2_GETH(is, id, h); }
#define IS2_ISEQ(is, id)  id = (id + 1) % IS2_N;
/* sequential picking */
#define IS2_PSEQ(is, id, h) { IS2_ISEQ(is, id); IS2_GETH(is, id, h); }

#define IS2_FLIP(is, id, h) { \
  is->M += (is->s[id] = -is->s[id]) * 2; \
  is->E += h * 2; }

#else

#define IS2_PICK(is, id, h)  id = is2_pick(is, &h)
#define IS2_FLIP(is, id, h)  is2_flip(is, id, h)
#endif /* IS2_LB */

#endif /* IS2_H__ */

