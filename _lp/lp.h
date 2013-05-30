#include "util.h"
#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef LP_H__
#define LP_H__
typedef struct {
  int n; /* number of non-basic variables */
  int m; /* number of constraints, i.e., basic variables */
  double *c;  /* c[1..n]: coefficients to objective function: -c[0]  + sum_i c[i] xi
               * array also the space holder for a[][] */
  double **a; /* a[1..m][0..n]: constraints */
  int *nid; /* indices for nonbasic variables */
  int *mid; /* indices for the basic variables */
} lpcore_t;

lpcore_t *lpcore_open(int n, int m);
void lpcore_close(lpcore_t *lp);
int lpcore_simplex(lpcore_t *lp);

#define lpcore_setc(lp, c) lpcore_seta(lp, c, 0)
INLINE void lpcore_seta(lpcore_t *lp, double *a, int j)
{
  int i;
  for (i = 0; i <= lp->n; i++) lp->a[j][i] = a[i];
  if (j == 0) lp->a[0][0] *= -1;
}

#define lpcore_printf(lp, fmt) lpcore_fprintf(lp, stdout, fmt)
INLINE void lpcore_fprintf(const lpcore_t *lp, FILE *fp, const char *fmt)
{
  int i, j;

  for (j = 0; j <= lp->m; j++) {
    if (j == 0) {
      fprintf(fp, "Objective function:\n");
    } else if (j == 1)
      fprintf(fp, "Coefficients:\n");
    for (i = 1; i <= lp->n; i++) {
      fprintf(fp, fmt, lp->a[j][i]);
      fprintf(fp, " x%-3d +", lp->nid[i]);
    }
    if (j == 0) {
      fprintf(fp, fmt, -lp->a[0][0]);
    } else {
      fprintf(fp, " x%-3d == ", lp->mid[j]);
      fprintf(fp, fmt, lp->a[j][0]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
}

#endif
