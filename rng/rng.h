#ifndef RNG_H__
#define RNG_H__

#define MTA_RAND  0
#define MTA_SAVE  1
#define MTA_DONE  2
#define MTA_FNAM  3
/* a random number in [0, 1) */
#define rnd0()        ((1.0/4294967296.0) * mtrand(MTA_RAND, 0, NULL))
#define mtsave(fn)    mtrand(MTA_SAVE, 0, fn)
#define mtfinish(fn)  mtrand(MTA_DONE, 0, fn)
#define mtsetfile(fn) mtrand(MTA_FNAM, 0, fn)
unsigned long mtrand(int action, unsigned long seed0, const char *fname);

double grand0(void);

#endif

