#include "util.h"
#ifndef HMC_H__
#define HMC_H__

typedef struct {
  int nd;
  real *x, *v, *f;
} hmc_t;


/* the macro versions may accept rv2_t and rv3_t */
#define hmc_push(h, x, v, f) \
  hmc_push_(h, (real *) x, (real *) v, (real *) f)

#define hmc_pop(h, x, v, f, vinv) \
  hmc_pop_(h, (real *) x, (real *) v, (real *) f, vinv)

#define hmc_open(nd, x, v, f) \
  hmc_open_(nd, (real *) x, (real *) v, (real *) f)

#define hmc_open2d(n, x, v, f) hmc_open(n*2, x, v, f)
#define hmc_open3d(n, x, v, f) hmc_open(n*3, x, v, f)

#define HMC_INVERTV 0x0001

/* save x, v, f as the new start point */
INLINE void hmc_push_(hmc_t *h, real *x, real *v, real *f)
{
  memcpy(h->x, x, h->nd * sizeof(real));
  memcpy(h->v, v, h->nd * sizeof(real));
  memcpy(h->f, f, h->nd * sizeof(real));
}


/* set x, v, f from the saved start point */
INLINE void hmc_pop_(hmc_t *h, real *x, real *v, real *f, unsigned flags)
{
  int i, nd = h->nd;

  /* invert velocities */
  if (flags & HMC_INVERTV) {
    for (i = 0; i < nd; i++) h->v[i] = -h->v[i];
  }
  memcpy(x, h->x, nd * sizeof(real));
  memcpy(v, h->v, nd * sizeof(real));
  memcpy(f, h->f, nd * sizeof(real));
}


/* open a hmc object from the start point */
hmc_t *hmc_open_(int nd, real *x, real *v, real *f)
{
  hmc_t *h;

  xnew(h, 1);
  h->nd = nd;
  xnew(h->x, nd);
  xnew(h->v, nd);
  xnew(h->f, nd);
  hmc_push(h, x, v, f);
  return h;
}

void hmc_close(hmc_t *h)
{
  free(h->x);
  free(h->v);
  free(h->f);
  free(h);
}

/* compute the reduce-momentum-flipping rate
 * J. Sohl-Dickstein, arXiv 1205.1939
 * the Metropolis probability of accepting x => x+ is min{1, exp(-dep)}
 * that of accepting x~ => x-~ is min{1, exp(-dem)}
 * where `~' means the momentum reversal operation */
int hmcgetredvflipr(double dem, double dep)
{
  die_if (dep <= 0, "the forward move must be rejected, dep %g\n", dep);
  /* 1. If the state x-, hence x-~, has higher energy than x+
   * dem >= dep, then x-~ is more likely rejected, and x~ should
   * inject probability into x; thus, there is no need to flip x
   * 2. We now have dem < dep.  We also have dep > 0, for x => x+
   * has been rejected. If dem <= 0, x~ => x-~ is definitely accepted,
   * then x should compensate the probability outflow of x~ and be flipped
   * 3. The only case left is 0 < dem < dep */
  if ( (dem >= dep)
    || (dem >= 0 && rnd0() >= (exp(dep-dem) - 1)/(exp(dep) - 1)) )
    return 0;
  return 1;
}

#endif

