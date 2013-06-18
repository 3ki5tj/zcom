#include "util.h"
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



lpcore_t *lpcore_open(int n, int m)
{
  lpcore_t *lp;
  int i, j;

  xnew(lp, 1);
  lp->n = n;
  lp->m = m;
  xnew(lp->c, (n+2)*(m+1));
  for (i = 0; i < (n+2)*(m+1); i++) lp->c[0] = 0.0f;
  xnew(lp->a, m+1);
  for (j = 0; j <= m; j++) /* map array indices */
    lp->a[j] = lp->c + j*(n + 2);
  xnew(lp->nid, n + 2);
  for (i = 0; i <= n + 1; i++) lp->nid[i] = i;
  xnew(lp->mid, m + 1);
  for (j = 0; j <= m; j++) lp->mid[j] = n + 1 + j;
  return lp;
}



void lpcore_close(lpcore_t *lp)
{
  free(lp->c);
  free(lp->a);
  free(lp->nid);
  free(lp->mid);
  free(lp);
}



/* mi: index of the entering (basic) variable, i.e. index in the objective function
 * mj: index of the leaving (basic) variable, i.e. index of equation */
static void lpcore_pivot(lpcore_t *lp, int mi, int mj)
{
  int i, j, m = lp->m, n = lp->n;
  double tmp, **a = lp->a;

  /* determine variable index */
  i = lp->mid[mj], lp->mid[mj] = lp->nid[mi], lp->nid[mi] = i;

  /* in row mj, solve mi in terms of mj and other nonbasic variables */
  for (tmp = 1.0/a[mj][mi], a[mj][mi] = 1.0, i = 0; i <= n; i++)
    a[mj][i] *= tmp;

  /* update the target function and other constraints */
  for (j = 0; j <= m; j++)
    if (j != mj)
      for (tmp = a[j][mi], a[j][mi] = 0.0, i = 0; i <= n; i++)
        a[j][i] -= tmp * a[mj][i];
}



/* solve linear programming */
static int lpcore_solve(lpcore_t *lp)
{
  int it, j, n = lp->n, m = lp->m;
  double **a = lp->a;
  double inc, tmp;
  int mi, mj;

  for (it = 0; ; it++) {
    /* 1. search the target function for a positive coefficient */
    for (mi = 1; mi <= n; mi++) if (a[0][mi] > 0) break;
    if (mi > n) break;

    /* 2. compute the most stringent constraint */
    for (mj = 0, inc = 1e30, j = 1; j <= m; j++)
      if (a[j][mi] > 0 && (tmp = a[j][0]/a[j][mi]) < inc) {
        inc = tmp;
        mj = j;
      }
    if (mj <= 0) return -1; /* unbounded */

    /* 3. swap non-basic variable mi with basic variable mj */
    lpcore_pivot(lp, mi, mj);
  }
  return 0;
}



/* update the objective function after reaching feasible solution */
static void lpcore_updobj(lpcore_t *lp, const double *c)
{
  int i, i1, j, idum, n = lp->n, m = lp->m;
  double **a = lp->a;

  /* clear the function first */
  for (a[0][0] = c[0], i = 1; i <= n; i++) a[0][i] = 0.;

  /* swap idum with the last nonbasic variable */
  idum = n + 1;
  for (i = 1; i <= n; i++)
    if (lp->nid[i] == idum) break;
  if (i != idum) {
    lp->nid[i] = lp->nid[idum];
    for (j = 1; j <= m; j++)
      a[j][i] = a[j][idum];
  }

  /* express each cj xj in terms of current nonbasic variables */
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= m; j++)
      if (lp->mid[j] == i) break;
    if (j <= m) {/* variable i is now basic, solve eq. j */
      for (i1 = 0; i1 <= n; i1++)
        a[0][i1] -= c[i]*a[j][i1];
    } else {
      for (i1 = 1; i1 <= n; i1++)
        if (lp->nid[i1] == i) break;
      a[0][i1] += c[i];
    }
  }
}



/* full simplex algorithm
 * return 0 on success, -1 if unbounded, -2 if no solution */
int lpcore_simplex(lpcore_t *lp)
{
  int i, idum, j, mj, m = lp->m, n = lp->n;
  double **a = lp->a, bmin, *c;

  for (bmin = 0, mj = 0, j = 1; j <= m; j++) {
    if (a[j][0] < bmin) {
      mj = j;
      bmin = a[j][0];
    }
  }
  if (mj == 0) return lpcore_solve(lp); /* there is a feasible solution */

  xnew(c, n + 1);
  for (i = 0; i <= n; i++) {
    c[i] = a[0][i];
    a[0][i] = 0;
  }
  idum = n + 1;
  for (j = 0; j <= m; j++) a[j][idum] = -1.0;
  lp->n = n + 1;
  /* search for a feasible solution */
  lpcore_pivot(lp, idum, mj);
  i = lpcore_solve(lp);
  if (i != 0) return i; /* unbounded */

  for (i = 1; i <= n+1; i++) /* see if the dummy variable is in the objective function */
    if (lp->nid[i] == idum) break;
  if (i <= n+1) { /* restore the original problem */
    lp->n = n;
    lpcore_updobj(lp, c);
    free(c);
    return lpcore_solve(lp);
  } else { /* no feasible solution */
    free(c);
    return -2;
  }
}
#endif


