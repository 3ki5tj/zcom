#ifndef IS_H__
#define IS_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
  int d, l, n;
  int E, M;
  int *s; /* 0 or 1 */
} is_t;

is_t   *is2_open(int l);
void    is2_close(is_t *is); 
int     is2_em(is_t *is);
int     is2_pick(is_t *is, int *nd);
int     is2_flip(is_t *is, int id, int nd);
int     is2_load(is_t *is, const char *fname);
int     is2_save(const is_t *is, const char *fname);
double  is2_exact(is_t *is, double beta, double *eav, double *cv);

#endif

