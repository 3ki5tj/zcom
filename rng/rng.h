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

#ifndef PRIu64
  #if defined(_MSC_VER) || defined(__BORLANDC__)
    #define PRIu32 "I32u"
    #define PRIu64 "I64u"
  #else
    #define PRIu32 "u"
    #define PRIu64 "llu"
  #endif
#endif

#define rand32()  mtrand()
#define rnd0()    ((1.0/4294967296.0) * rand32()) /* double, [0, 1) */

#define MTFILE    "MTSEED"  /* default file */
#define MTSEED    5489UL    /* default seed */
int mtsave(const char *fname);
int mtload(const char *fname, uint32_t seed);
uint32_t mtrand(void);
double grand0(void);

#endif
