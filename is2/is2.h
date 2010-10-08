#include "rng.h"
#ifndef IS2_H__
#define IS2_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
  int d, l, n;
  int E, M;
  int *s; /* 0 or 1 */
} is_t;

int     is2_em(is_t *is);
int     is2_load(is_t *is, const char *fname);
int     is2_save(const is_t *is, const char *fname);
double  is2_exact(is_t *is, double beta, double *eav, double *cv);
int     is2_pick(const is_t *is, int *nd);
int     is2_flip(is_t *is, int id, int nd);

#ifdef  IS2_L
#define L_  IS2_L
#define N_  (L_ * L_)

static int S_[N_];

#define IS2_PICK(is, id, nd) { \
  int i, j; \
  id = (int)(rnd0() * N_); \
  i = id / L_, j = id % L_; \
  nd =  S_[i * L_ + (j+L_-1)%L_] + S_[i * L_ + (j+1)%L_] \
      + S_[(i-1+L_)%L_ * L_ + j] + S_[(i+1)%L_ * L_ + j] -2 ; \
  if (!S_[id]) nd = -nd; }

#define IS2_FLIP(is, id, nd) \
  S_[id] = !S_[id]; \
  is->E += nd*4;

  //is->M += (S_[id] = !S_[id]) ? 2 : -2; \

#endif /* IS2_L */

/* initialize an lxl Ising model */
static __inline is_t *is2_open(int l)
{
  int i, n;
  is_t *is;

  if ((is = calloc(1, sizeof(*is))) == NULL){
    fprintf(stderr, "no memory for is.\n");
    return NULL;
  }
  is->d = 2;
  is->l = l;
  is->n = n = l*l;
#ifdef IS2_L
  assert(l == IS2_L);
  is->s = S_;
#else
  if ((is->s = malloc(sizeof(is->s[0])*n)) == NULL) {
    fprintf(stderr, "no memory for spin, %dx%d\n", l, l);
    return NULL;
  }
#endif
  for (i = 0; i < n; i++) is->s[i] = 0;
  is->M = -n;
  is->E = -2*n;
  return is;
}

static __inline void is2_close(is_t *is) 
{
  if (is != NULL) {
#ifndef IS2_L
    free(is->s);
#endif
    free(is);
  }
}

#endif /* IS2_H__ */

