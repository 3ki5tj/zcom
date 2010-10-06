/* stuff before the first #ifndef is not imported to zcom.h */
#ifndef RNG_C__
#define RNG_C__
#include "rng.h"

/* Mersenne Twister was developped by Makoto Matsumoto and Takuji Nishimura */
#define N_MT 624
#define M_MT 397
#define UMASK_MT 0x80000000UL /* most significant w-r bits */
#define LMASK_MT 0x7fffffffUL /* least significant r bits */

/* save the current mt state to file */
static int mtsave_(int idx, unsigned long arr[], const char *fname)
{
  FILE *fp;
  int k;

  if (idx < 0) return 1; /* RNG was never used, so it cannot be saved */
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot save to %s.\n", fname);
    return 1;
  }
  fprintf(fp, "MTSEED\n%d\n", idx);
  for (k = 0; k < N_MT; k++) fprintf(fp, "%lu\n", arr[k]);
  fclose(fp);
  return 0;
}

/* load mt state from `fname', or if it fails, use `seed' to initialize mt  */
static void mtload_(int *pindex, unsigned long arr[], const char *fname,
    unsigned long seed)
{
  char s[32];
  int k, z, err = 1;
  FILE *fp;

  if ((fp = fopen(fname, "r")) != NULL) { /* try to load from file */
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "%s is empty\n", fname);
    } else if (strncmp(s, "MTSEED", 6) != 0) { /* to check the first line */
      fprintf(stderr, "mtrand: corrupted file.\n");
    } else if (fscanf(fp, "%d", pindex) != 1) {
      fprintf(stderr, "no index in %s\n", fname);
    } else {
      if (*pindex < 0) *pindex = N_MT; /* request updating */
      for (z = 1, k = 0; k < N_MT; k++) {
        if (fscanf(fp, "%lu", &arr[k]) != 1) break;
        if (arr[k] != 0) z = 0; /* a non-zero number */
      }
      if (k != N_MT) fprintf(stderr, "%s incomplete %d/%d\n", fname, k, N_MT);
      else err = z; /* clear error, if array is nonzero */
    }
    fclose(fp);
  }

  if (err) { /* initialize from seed */
    if (seed == 0) seed = MTSEED;
    arr[0] = seed & 0xffffffffUL;
    for (k = 1; k < N_MT; k++) /* the final mask is for 64-bit machines */
      arr[k] = (1812433253UL * (arr[k-1] ^ (arr[k-1]>>30)) + k) & 0xffffffffUL;
    *pindex = N_MT; /* request updating */
  }
}

/* return an unsigned random number */
static unsigned long  mtrand_(int *pindex, unsigned long arr[])
{
  const unsigned long mag01[2] = { 0, 0x9908b0dfUL }; /* MATRIX_A */
  unsigned long x;
  int k, kp;

  if (*pindex >= N_MT) { /* generate N_MT words at one time */
    for (k = 0; k < N_MT - M_MT; k++) {
      x = (arr[k] & UMASK_MT) | (arr[k+1] & LMASK_MT);
      arr[k] = arr[k+M_MT] ^ (x>>1) ^ mag01[x&1UL];
    }
    for (; k < N_MT; k++) {
      if ((kp=k+1) == N_MT) kp = 0;
      x = (arr[k] & UMASK_MT) | (arr[kp] & LMASK_MT);
      arr[k] = arr[k+(M_MT-N_MT)] ^ (x>>1) ^ mag01[x&1UL];
    }
    *pindex = 0;
  }
  x = arr[ (*pindex)++ ];
  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680UL;
  x ^= (x << 15) & 0xefc60000UL;
  x ^= (x >> 18);
  return x;
}

/* random number generator, `cmd' can be:
 * MTA_RAND: generate a 32-bit random number, try to load state from file
 *           at the beginning, and if unsuccessful, use `seed' to init.
 * MTA_SAVE: save the current state to file */
unsigned long mtrand(int cmd, unsigned long seed, const char *fnm)
{
  static unsigned long arr[N_MT]; /* array for the mt state vector */
  static int idx = -1; /* index in arr, -1: uninitialized */
  const char *fname = (fnm != NULL) ? fnm : MTFILE;

  switch (cmd) {
  case MTA_RAND:   /* generate random number */
    if (idx < 0) mtload_(&idx, arr, fname, seed);
    return mtrand_(&idx, arr);
  case MTA_SAVE:   /* save the current state  */
    return mtsave_(idx, arr, fname);
  default:
    fprintf(stderr, "mtrand: unknown cmd %d\n", cmd);
  }
  return 0;
}
#undef N_MT
#undef M_MT
#undef UMASK_MT
#undef LMASK_MT

/* Gaussian distribution with zero mean and unit variance */
double grand0(void)
{
  const double tol = 1e-15; /* 1e-13 allows fac of 1.0/(6e-8) */
  static double save = 0.0;
  static int stock = 0;
  double fac, r2, x, y;

  if (stock) { stock = 0; return save; }
  for (;;) {
    x = 2.0 * rnd0()-1.0;
    y = 2.0 * rnd0()-1.0;
    r2 = x * x + y * y;
    /* to make sure it's inside a unit circle, and
       r2 can be a denominator */
    if (r2 < 1.0 && r2 > tol) break;
  }
  fac   = sqrt(-2.0*log(r2)/r2);
  stock = 1;
  save  = y * fac;
  return x * fac;
}

#endif
