#define ZCINLINE __inline static
#ifndef AV_H__
#define AV_H__
#include <stdio.h>
#include <math.h>

typedef struct {
  double s, sx, sx2; /* sum, sum x, variance */
} av_t;

void av_add(av_t *av, double w, double x);
ZCINLINE void av_clear(av_t *av) { av->s = av->sx = av->sx2 = 0; }
ZCINLINE double av_getave(const av_t *av) { return av->s > 0 ? av->sx/av->s : 0; }
ZCINLINE double av_getvar(const av_t *av) { return av->s > 0 ? av->sx2/av->s : 0; }
ZCINLINE double av_getdev(const av_t *av) { return av->s > 0 ? sqrt(av_getvar(av)) : 0; }

/* add a new value to av_t with a weight `w' */
ZCINLINE void av_addw(av_t *av, double x, double w)
{
  double s, sx;

  av->s = (s = av->s) + w;
  av->sx = (sx = av->sx) + x*w;
  if (s <= 0.0) return;
  av->sx2 += (x - sx/s)*(x - av->sx/av->s)*w;
}
#define av_add(av, x) av_addw(av, x, 1)


#endif
