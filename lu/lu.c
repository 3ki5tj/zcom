#ifndef LU_C__
#define LU_C__

/* LU decomposition part  */
#include "lu.h"

/* solve A x = b by L U decomposition */
int lusolve(real *a, real *b, int n)
{
  int i, j, k, imax = 0;
  real x, max;
  const real mintol = 1e-16; /* absolute minimal value for a pivot */

  for (i = 0; i < n; i++) {  /* normalize each equation */
    for (max = 0.0, j = 0; j < n; j++)
      if ((x = fabs(a[i*n+j])) > max) 
        max = x;
    if (max < mintol) {
      return 1;
    }
    for (x = 1.0/max, j = 0; j < n; j++) 
      a[i*n+j] *= x;
    b[i] *= x;
  }

  /* step 1: A = L U, column by column */
  for (j = 0; j < n; j++) {
    /* matrix U */
    for (i = 0; i < j; i++) {
      for (x = a[i*n+j], k = 0; k < i; k++)
        x -= a[i*n+k]*a[k*n+j];
      a[i*n+j] = x;
    }

    /* matrix L, diagonal of L are 1 */
    max = 0.0;
    for (i = j; i < n; i++) {
      for (x = a[i*n+j], k = 0; k < j; k++) 
        x -= a[i*n+k]*a[k*n+j];
      a[i*n+j] = x;
      if (fabs(x) >= max) { 
        max = fabs(x);
        imax = i;
      }
    }

    if (j != imax) { /* swap the pivot row with the jth row */
      for (k = 0; k < n; k++)
        x = a[imax*n+k], a[imax*n+k] = a[j*n+k], a[j*n+k] = x;
      x = b[imax]; b[imax] = b[j]; b[j] = x;
    }
    if (fabs(a[j*n+j]) < mintol)
      return 2;
    /* divide by the pivot element, for the L matrix */
    if (j != n-1)
      for (x = 1.0/a[j*n+j], i = j+1; i < n; i++)
        a[i*n+j] *= x;
  }

  /* step2: solve the equation L U x = b */
  for (i = 0; i < n; i++) { /* L y = b */
    x = b[i];
    for (j = 0; j < i; j++) x -= a[i*n+j]*b[j];
    b[i] = x;
  }
  for (i = n-1; i >= 0; i--) { /* U x = y. */
    x = b[i];
    for (j = i+1; j < n; j++) x -= a[i*n+j]*b[j];
    b[i] = x/a[i*n+i];
  }
  return 0;
}

#endif

