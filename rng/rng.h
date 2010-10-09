#ifndef RNG_H__
#define RNG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef UINT32  /* 32-bit unsigned int */
#define UINT32  unsigned
#define UI32FMT "%u"
#endif

#define rand32()  mtrand()
#define rnd0()    ((1.0/4294967296.0) * rand32()) /* double, [0, 1) */

#define MTFILE    "MTSEED"  /* default file */
#define MTSEED    5489UL    /* default seed */
int mtsave(const char *fname);
int mtload(const char *fname, UINT32 seed);
UINT32 mtrand(void);
double grand0(void);

#endif
