#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef RNG_H__
#define RNG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
  #include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
  typedef unsigned uint32_t;
  typedef unsigned __int64 uint64_t;
#else
  #include <inttypes.h>
#endif

#ifndef PRIu32
  #if (defined(_MSC_VER) && (_MSC_VER >= 1300)) || defined(__BORLANDC__)
    #define PRIu32 "I32u"
  #else
    #define PRIu32 "u"
  #endif
#endif

#ifndef PRIu64
  #if defined(_MSC_VER) || defined(__BORLANDC__)
    #define PRIu64 "I64u"
  #else
    #define PRIu64 "llu"
  #endif
#endif

#define rand32()  mtrand()
#define rnd0()    ((1.0/4294967296.0) * rand32()) /* double, [0, 1) */

#define MTFILE    "MTSEED"  /* default file */
#define MTSEED    5489UL    /* default seed */
INLINE int mtsave(const char *fname);
INLINE int mtload(const char *fname, uint32_t seed);
INLINE uint32_t mtrand(void);
INLINE double grand0(void);

/* metropolis acceptance probability rnd0() < exp(r), assuming r > 0 */
INLINE int metroacc0(double r) { r = exp(r); return rnd0() < r; }

/* metropolis acceptance probability rnd0() < exp(- bet * de), assuming bet > 0 
 * defined as a macro, in case r is an integer */
#define metroacc1(de, bet) ((de <= 0) ? 1 : metroacc0(- bet * de))

#endif
