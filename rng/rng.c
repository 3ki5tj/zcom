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
INLINE uint32_t mtrand(void)
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
INLINE double grand0(void)
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


/* return a random number that satisfies a gamma distribution
   p(x) = x^(k - 1) e^(-x) / (k - 1)!. */
INLINE double randgam(int k)
{
  int i;
  double x, k1 = k - 1, r, y, v1, v2, s;

  if (k < 0) { printf("randgam: k %d must be positive\n", k); return 0.; }
  if (k == 0) return 0.; /* nothing */
  if (k <= 7) { /* adding numbers of exponential distribution */
    /* exp(- x1 - x2 - x3 - x4) dx1 dx2 dx3 dx4 */
    for (x = 1.0, i = 0; i < k; i++)
      x *= rnd0();
    return -log(x);
  }

  /* generate gamma distribution by the rejection method */
  for(;;) {
    /* generate lorentz distribution, centered at k1, width is sqrt(2.0*k - 1)
     p(y) = 1/pi/(1 + y^2), x = y*w + k1, w = sqrt(2.0*k - 1) */
    for (;;) { /* get a unit circle */
      v1 = 2.0*rnd0() - 1.0;
      v2 = 2.0*rnd0() - 1.0;
      if (v1*v1 + v2*v2 <= 1.0) {
        y = v2/v1; /* tan */
        s = sqrt(2.0*k - 1);
        x = s*y + k1;
        if (x > 0.0) break; /* drop the negative value */
      }
    }
    /* compare with the gamma distribution
       r peaks at x = k1, where, y = 0 and r = 1 */
    r = (1.0 + y*y)*exp(k1*log(x/k1) - x + k1);
    if (rnd0() <= r) break;
  }

  return x;
}

/* return the sum of the square of Gaussian random numbers  */
INLINE double randgausssum(int n)
{
  double x, r; 
  if (n <= 0) return 0.0;
  x = 2.0*randgam(n/2);
  if (n % 2) { r = grand0(); x += r*r; }
  return x;
}

/* return randomly oriented vector on a uniform sphere */
INLINE double *randor3d(double *v)
{
  double a, b, sq, s;

  do { /* projection on the x-y plane */
    a = 2 * rnd0() - 1;
    b = 2 * rnd0() - 1;
    sq = a * a + b * b;
  } while (sq > 1);

  s = 2.0 * sqrt(1 - sq);

  /* x^2 + y^2 = 4 [1 - (a^2 + b^2)] (a^2 + b^2)
   * z^2 = 1 - 2 (a^2 + b^2), so x^2 + y^2 + z^2 = 1 */
  v[0] = a * s;
  v[1] = b * s;
  v[2] = 1 - 2*sq;

  return v;
}

#endif
