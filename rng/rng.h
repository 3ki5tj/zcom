#ifndef RNG_H__
#define RNG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>

#define MTFILE "MTSEED" /* default file */
#define MTSEED 5489UL /* default seed */
#define rnd0() ((1.0/4294967296.0) * mtrand()) /* double, [0, 1) */

int mtsave(const char *fname);
int mtload(const char *fname, unsigned long seed);
unsigned long mtrand(void);
double grand0(void);

#endif
