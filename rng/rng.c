#include "ss.h"

/* stuff before the first #ifndef is not imported to zcom.h */
#ifndef RNG_C__
#define RNG_C__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rng.h"

/* The following random number generator algorithm Mersenne Twister
   was developped by Makoto Matsumoto and Takuji Nishimura */

/* Period parameters */
#define N_MT 624
#define M_MT 397
#define UMASK_MT 0x80000000UL /* most significant w-r bits */
#define LMASK_MT 0x7fffffffUL /* least significant r bits */

/* save the current mt state to file */
static int mtsavestate_(int idx, unsigned long arr[], const char *fname)
{
  FILE *fp;
  int k;

  if (idx < 0) /* RNG was never used, so it cannot be saved */
    return 1;
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "mtrand: cannot save state to %s.\n", fname);
    return 1;
  }
  fprintf(fp, "MTSEED\n%d\n", idx);
  for (k = 0; k < N_MT; k++)
    fprintf(fp, "%lu\n", arr[k]);
  fclose(fp);
  return 0;
}

/* read mt state from file, or if it fails, use seed0 to initialize mt  */
static void mtloadstate_(int *pindex, unsigned long arr[], const char *fname,
    unsigned long seed0)
{
  unsigned long mtseed0 = 5489UL;
  char *s = NULL, err = 1;
  int nz, k;
  FILE *fp;

  /* here we try to initialize the array from file */
  if ((fp = fopen(fname, "r")) != NULL) {
    if (ssfgets(s, NULL, fp) == NULL) {
      fprintf(stderr, "mtrand: cannot read the first line of %s, "
          "probably an empty file.\n", fname);
      goto CLOSEFILE;
    }

    if (strncmp(s, "MTSEED", 6) != 0) { /* to check the first line */
      fprintf(stderr, "mtrand: corrupted seed file.\n");
    } else {
      if (fscanf(fp, "%d", pindex) != 1) {
        fprintf(stderr, "mtrand: no index found.\n");
        goto CLOSEFILE;
      }
      if (*pindex < 0) /* request updating */
        *pindex = N_MT;
      for (k = 0; k < N_MT; k++) {
        if (fscanf(fp, "%lu", &arr[k]) != 1)
          break;
        if (arr[k] != 0) /* a non-zero number */
          nz = 1;
      }
      if (k != N_MT)
        fprintf(stderr, "mtrand: seed file is incomplete, got %d numbers.\n", k);
      else
        err = 0; /* clear error */
    }
CLOSEFILE:
    fclose(fp);
    ssdel(s);
  }

  if (!nz) /* last item read is zero */
    err = 1;

  if (err) {
    if (seed0 != 0)
      mtseed0 = seed0;
    arr[0] = mtseed0 & 0xffffffffUL;
    for (k = 1; k < N_MT; k++) {
      arr[k] = (1812433253UL * (arr[k-1] ^ (arr[k-1]>>30)) + k) & 0xffffffffUL;
    } /* the masking step is for 64-bit machines */
    *pindex = N_MT; /* request updating */
  }
}

/* return an unsigned random number */
static unsigned long  mtrandlow_(int *pindex, unsigned long arr[])
{
  const unsigned long mag01[2] = { 0, 0x9908b0dfUL }; /* MATRIX_A */
  unsigned long x;
  int k;

  if (*pindex >= N_MT) { /* generate N_MT words at one time */
    for (k = 0; k < N_MT - M_MT; k++) {
      x = (arr[k] & UMASK_MT) | (arr[k+1] & LMASK_MT);
      arr[k] = arr[k+M_MT] ^ (x>>1) ^ mag01[x&1UL];
    }
    for (; k < N_MT - 1; k++) {
      x = (arr[k] & UMASK_MT) | (arr[k+1] & LMASK_MT);
      arr[k] = arr[k+(M_MT-N_MT)] ^ (x>>1) ^ mag01[x&1UL];
    }
    x = (arr[N_MT-1]&UMASK_MT) | (arr[0]&LMASK_MT);
    arr[N_MT-1] = arr[M_MT-1] ^ (x>>1) ^ mag01[x&1UL];
    *pindex = 0;
  }

  x = arr[ (*pindex)++ ];

  /* Tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680UL;
  x ^= (x << 15) & 0xefc60000UL;
  x ^= (x >> 18);

  return x;
}

/* The first argument defines a command
 * MTA_RAND: to generate random number unsigned long
 *           for the first call, we try to load the state from file `fname',
 *           if it is unsuccessful, seed0, if unzero, is used to initialize the RNG
 * MTA_SAVE: save the current state to `fname'
 * MTA_DONE: save state and finish
 * MTA_FNAM: change the default file name
 * */
unsigned long mtrand(int action, unsigned long seed0, const char *fname0)
{
  static unsigned long arr[N_MT]; /* array for the mt state vector */
  static int idx = -1; /* index in arr, -1: uninitialized */
  static char *fname = NULL;

  if (fname0 != NULL) /* if fname is given (not NULL), always use it */
    sscpy(fname, fname0);
  else if (fname == NULL) /* otherwise,  use the default */
    sscpy(fname, "MTSEED");

  switch (action) {
  case MTA_RAND:   /* generate random number */
    if (idx < 0)
      mtloadstate_(&idx, arr, fname, seed0);
    return mtrandlow_(&idx, arr);

  case MTA_SAVE:   /* save the current state  */
    mtsavestate_(idx, arr, fname);
    break;

  case MTA_DONE:   /* finish up */
    mtsavestate_(idx, arr, fname);
    ssdelete(fname);
    break;

  case MTA_FNAM:   /* set fname, done */
    break;

  default:
    fprintf(stderr, "mtrand: unknown action %d\n", action);
  }
  return 0;
}

#undef N_MT
#undef M_MT
#undef UMASK_MT
#undef LMASK_MT


/* Gaussian distribution with zero mean and unit variance,
   Algorithm adapted from Numerical recipes. */
double grand0(void)
{
  const double  tol = 1e-15; /* 1e-13 allows fac of 1.0/(6e-8) */
  static double save = 0.0;
  static int    stock = 0;
  double fac, r2, x, y;

  if (stock) {
    stock = 0;
    return save;
  }

  for (;;) {
    x = 2.0 * rnd0()-1.0;
    y = 2.0 * rnd0()-1.0;
    r2 = x * x + y * y;
    /* to make sure it's inside a unit circle, and
       r2 can be a denominator */
    if (r2 < 1.0 && r2 > tol)
      break;
  }
  fac   = sqrt(-2.0*log(r2)/r2);
  stock = 1;
  save  = y * fac;
  return x * fac;
}

#endif

