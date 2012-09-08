#ifndef RNG_C__
#define RNG_C__
#include "rng.h"

/* Mersenne Twister was developped by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000UL /* most significant w-r bits */
#define MT_LMASK 0x7fffffffUL /* least significant r bits */

int mtidx_ = -1; /* index in mt_, -1: uninitialized */
uint32_t mt_[MT_N]; /* array for the mt state vector */

/* save the current mt state to file */
INLINE int mtsave(const char *fname)
{
  FILE *fp;
  int k;

  if (mtidx_ < 0) return 1; /* RNG was never used, so it cannot be saved */
  if (fname == NULL) fname = MTFILE;
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot save to %s.\n", fname);
    return 1;
  }
  fprintf(fp, "MTSEED\n%d\n", mtidx_);
  for (k = 0; k < MT_N; k++) fprintf(fp, "%"PRIu32"\n", mt_[k]);
  fclose(fp);
  return 0;
}

/* load mt state from `fname', or if it fails, use `seed' to initialize mt  */
INLINE int mtload(const char *fname, uint32_t seed)
{
  static char s[64];
  int k, z, err = 1;
  FILE *fp;
  
  if (fname == NULL) fname = MTFILE;
  if ((fp = fopen(fname, "r")) != NULL) { /* try to load from file */
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "%s is empty\n", fname);
    } else if (strncmp(s, "MTSEED", 6) != 0) { /* to check the first line */
      fprintf(stderr, "mtrand: corrupted file.\n");
    } else if (fscanf(fp, "%d", &mtidx_) != 1) {
      fprintf(stderr, "no index in %s\n", fname);
    } else {
      if (mtidx_ < 0) mtidx_ = MT_N; /* request updating */
      for (z = 1, k = 0; k < MT_N; k++) {
        if (fscanf(fp, "%"PRIu32, &mt_[k]) != 1) break;
        if (mt_[k] != 0) z = 0; /* a non-zero number */
      }
      if (k != MT_N) fprintf(stderr, "%s incomplete %d/%d\n", fname, k, MT_N);
      else err = z; /* clear error, if array is nonzero */
    }
    fclose(fp);
  }

  if (err) { /* initialize from seed */
    if (seed == 0) seed = MTSEED;
    mt_[0] = seed & 0xffffffffUL;
    for (k = 1; k < MT_N; k++) /* the final mask is for 64-bit machines */
      mt_[k] = (1812433253UL * (mt_[k-1] ^ (mt_[k-1]>>30)) + k) & 0xffffffffUL;
    mtidx_ = MT_N; /* request updating */
  }
  return (mtidx_ < 0);
}

/* return an unsigned random number */
uint32_t mtrand(void)
{
  uint32_t x;
  static const uint32_t mag01[2] = {0, 0x9908b0dfUL}; /* MATRIX_A */
  int k;

  if (mtidx_ < 0) mtload(NULL, 0);
  if (mtidx_ >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mt_[k] & MT_UMASK) | (mt_[k+1] & MT_LMASK);
      mt_[k] = mt_[k+MT_M] ^ (x>>1) ^ mag01[x&1UL];
    }
    for (; k < MT_N-1; k++) {
      x = (mt_[k] & MT_UMASK) | (mt_[k+1] & MT_LMASK);
      mt_[k] = mt_[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1UL];
    }
    x = (mt_[MT_N-1] & MT_UMASK) | (mt_[0] & MT_LMASK);
    mt_[MT_N-1] = mt_[MT_M-1] ^ (x>>1) ^ mag01[x&1UL];
    mtidx_ = 0;
  }
  x = mt_[ mtidx_++ ];
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
double grand0(void)
{
  double x, y, u, v, q;
  do {
    u = 1 - rnd0();
    v = 1.7156*(rnd0() - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}

#endif
