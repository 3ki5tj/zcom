#include "util.h"
#ifndef RNG_H__
#define RNG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>

#define MTFILE    "MTSEED"  /* default file */
#define MTSEED    5489UL    /* default seed */

/* Mersenne Twister was developed by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000UL /* most significant w-r bits */
#define MT_LMASK 0x7fffffffUL /* least significant r bits */

typedef struct {
  int idx;
  uint32_t arr[MT_N];
} mtrng_t;


/* local copy for convenience */
mtrng_t mrstock_ = {-1, {0}}; /* index in mr_, -1: uninitialized */

mtrng_t *mr_ = &mrstock_;
/* by default `mr_' is pointed to `mtstock_', but we replace it
 * by a thread private version by calling mtmkprivate(seed)
 * this trick allows the different random numbers generated
 * from different threads even by calling the default versions
 * of the functions, e.g., rand01() */
#ifdef _OPENMP
#pragma omp threadprivate(mr_)
#endif



/* in the default scrambling, we create a new RNG to be safe */
#define mtscramble(seed) { \
  if (mr_ == &mrstock_) mr_ = mtrng_open0(); \
  mtrng_scramble(mr_, seed); }

/* free the default RNG */
#define mtclosedef() { if (mr_ != &mrstock_) mtrng_close(mr_); }


/* default versions */
#define mtsave(fn)          mtrng_save(mr_, fn)
#define mtload(fn, seed)    mtrng_load(mr_, fn, seed)
#define mtrand()            mtrng_rand(mr_)
#define metroacc0(r)        mtrng_metroacc0(mr_, r)
#define metroacc1(de, bet)  mtrng_metroacc1(mr_, de, bet)

/* old aliases */
#define rnd0()              rand01()
#define rnd(a, b)           randunif(a, b)
#define grand0()            randgaus()



INLINE mtrng_t *mtrng_open0(void)
{
  mtrng_t *mr;

  xnew(mr, 1);
  mr->idx = -1;
  return mr;
}



INLINE void mtrng_close(mtrng_t *mr)
{
  free(mr);
}



/* save the current state to file */
INLINE int mtrng_save(mtrng_t *mr, const char *fn)
{
  FILE *fp;
  int k;

  if (mr->idx < 0) return 1; /* RNG was never used, so it cannot be saved */
  if (fn == NULL) fn = MTFILE;
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "MTSEED\n%d\n", mr->idx);
  for (k = 0; k < MT_N; k++)
    fprintf(fp, "%" PRIu32 "\n", mr->arr[k]);
  fclose(fp);
  return 0;
}


/* randomize the array `mr_' */
INLINE void mtrng_scramble(mtrng_t *mr, uint32_t seed)
{
  int k;

  mr->arr[0] = ((seed + MTSEED) * 314159265ul + 271828183ul) & 0xfffffffful;
  for (k = 1; k < MT_N; k++) { /* the final mask is for 64-bit machines */
    mr->arr[k] = 1812433253ul * (mr->arr[k - 1] ^ (mr->arr[k - 1] >> 30)) + k;
    /* mr->arr[k] = (mr->arr[k] + seed) * 22695477ul + 1ul; */
    mr->arr[k] = ((mr->arr[k] + seed) * 314159265ul + 1ul) & 0xfffffffful;
  }
  mr->idx = MT_N; /* request for an update */
}



INLINE mtrng_t *mtrng_open(uint32_t seed)
{
  mtrng_t *mr = mtrng_open0();
  mtrng_scramble(mr, seed);
  return mr;
}



/* load mr state from `fn', or if it fails, use `seed' to initialize mr  */
INLINE int mtrng_load(mtrng_t *mr, const char *fn, uint32_t seed)
{
  char s[64];
  int k, z, err = 1;
  FILE *fp;

  if (fn == NULL) fn = MTFILE;
  if ((fp = fopen(fn, "r")) != NULL) { /* try to load from file */
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "%s is empty\n", fn);
    } else if (strncmp(s, "MTSEED", 6) != 0) { /* to check the first line */
      fprintf(stderr, "%s corrupted\n", fn);
    } else if (fscanf(fp, "%d", &mr->idx) != 1) {
      fprintf(stderr, "no index in %s\n", fn);
    } else {
      if (mr->idx < 0) mr->idx = MT_N; /* request updating */
      for (z = 1, k = 0; k < MT_N; k++) {
        if (fscanf(fp, "%" SCNu32, &mr->arr[k]) != 1) break;
        if (mr->arr[k] != 0) z = 0; /* a non-zero number */
      }
      if (k != MT_N) fprintf(stderr, "%s incomplete %d/%d\n", fn, k, MT_N);
      else err = z; /* clear error, if array is nonzero */
    }
    fclose(fp);
  }

  if (err) mtrng_scramble(mr, seed);
  return (mr->idx < 0);
}



#define mtrng_rand32(mr) mtrng_rand(mr)
/* must be double to avoid the round-off error that gives >= 1 result */
#define mtrng_rand01(mr) (mtrng_rand32(mr) * (1./4294967296.0)) /* double, [0, 1) */
#define mtrng_randunif(mr, a, b) ((a) + ((b) - (a)) * mtrng_rand01(mr)) /* double, [a, b) */

/* return an unsigned random number */
INLINE uint32_t mtrng_rand(mtrng_t *mr)
{
  static const uint32_t mag01[2] = {0, 0x9908b0dfUL}; /* MATRIX_A */
#ifdef _OPENMP
#pragma omp threadprivate(mag01)
#endif
  uint32_t x;
  int k;

  if (mr->idx < 0) mtrng_load(mr, NULL, 0);
  if (mr->idx >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mr->arr[k] & MT_UMASK) | (mr->arr[k+1] & MT_LMASK);
      mr->arr[k] = mr->arr[k+MT_M] ^ (x>>1) ^ mag01[x&1UL];
    }
    for (; k < MT_N-1; k++) {
      x = (mr->arr[k] & MT_UMASK) | (mr->arr[k+1] & MT_LMASK);
      mr->arr[k] = mr->arr[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1UL];
    }
    x = (mr->arr[MT_N-1] & MT_UMASK) | (mr->arr[0] & MT_LMASK);
    mr->arr[MT_N-1] = mr->arr[MT_M-1] ^ (x>>1) ^ mag01[x&1UL];
    mr->idx = 0;
  }
  x = mr->arr[ mr->idx++ ];
  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680UL;
  x ^= (x << 15) & 0xefc60000UL;
  x ^= (x >> 18);
  return x;
}



#undef MT_N
#undef MT_M
#undef MT_UMASK
#undef MT_LMASK



/* Gaussian distribution with zero mean and unit variance
 * using ratio method */
INLINE double mtrng_randgaus(mtrng_t *mr)
{
  double x, y, u, v, q;
  do {
    u = 1 - mtrng_rand01(mr);
    v = 1.7156*(mtrng_rand01(mr) - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}



/* return a random number that satisfies the gamma distribution
 * p(x) = x^(k - 1) exp(-x) / Gamma(k) */
INLINE double mtrng_randgam(mtrng_t *mr, double k)
{
  int lt1 = 0;
  double a, b, x, v, u;

  if (k <= 0) {
    fprintf(stderr, "mtrng_randgam: k %g must be positive\n", k);
    return 0.;
  }
  if ( k < 1 ) {
    lt1 = 1;
    k += 1;
  }
  a = k - 1./3;
  b = 1./3/sqrt(a);

  for ( ; ; ) {
    do {
      x = mtrng_randgaus(mr);
      v = 1 + b * x;
    } while ( v <= 0 );
    v *= v * v;
    x *= x;
    u = mtrng_rand01(mr);
    if ( u <= 1 - 0.331 * x * x ) break;
    u = log(u);
    if ( u <= 0.5 * x + a * (1 - v + log(v)) ) break;
  }

  x = a * v;
  if ( lt1 ) x *= pow(1 - mtrng_rand01(mr), 1./(k - 1));
  return x;
}



/* return a random number that satisfies the chi-squared distribution,
 * which is the sum of the squares n Gaussian random numbers */
INLINE double mtrng_randchisqr(mtrng_t *mr, double n)
{
  return 2 * mtrng_randgam(mr, n*.5);
}



/* random pair index (i, j) */
INLINE int mtrng_randpair(mtrng_t *mr, int n, int *j)
{
  int pid = (int) (mtrng_rand01(mr) * n * (n - 1)), i;
  i = pid / (n - 1);
  *j = pid - i * (n - 1);
  if (*j >= i) (*j)++;
  return i;
}



/* the following are declared as functions instead of macros
 * so they can be used as function pointers */

INLINE uint32_t rand32(void)
{
  return mtrng_rand32(mr_);
}

INLINE double rand01(void)
{
  return mtrng_rand01(mr_);
}

INLINE double randunif(double a, double b)
{
  return mtrng_randunif(mr_, a, b);
}

INLINE double randgaus(void)
{
  return mtrng_randgaus(mr_);
}

INLINE double randgam(double k)
{
  return mtrng_randgam(mr_, k);
}

INLINE double randchisqr(double n)
{
  return mtrng_randchisqr(mr_, n);
}

INLINE int randpair(int n, int *j)
{
  return mtrng_randpair(mr_, n, j);
}



/* Metropolis acceptance probability rand01() < exp(- bet * de), assuming bet > 0
 * defined as a macro, in case r is an integer */
#define mtrng_metroacc1(mr, de, bet) \
  ((de <= 0) ? 1 : mtrng_metroacc0(mr, -bet * de))

/* Metropolis acceptance probability rand01() < exp(r), assuming r > 0 */
INLINE int mtrng_metroacc0(mtrng_t *mr, double r)
{
  r = exp(r);
  return mtrng_rand01(mr) < r;
}



#endif

