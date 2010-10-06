#ifndef RNG_H__
#define RNG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>
#define MTFILE "MTSEED" /* default file */
#define MTSEED 5489UL /* default seed */
#define MTA_RAND  0
#define MTA_SAVE  1
#define rand32()      mtrand(MTA_RAND, 0, NULL)
#define rnd0()        ((1.0/4294967296.0) * rand32()) /* double, [0, 1) */
#define mtsave(fn)    mtrand(MTA_SAVE, 0, fn)
unsigned long mtrand(int action, unsigned long seed0, const char *fname);
double grand0(void);

#endif

