#include <stdio.h>

#ifdef REF
#include "zcomref.h"
#else
#include "rng.h"
#endif

int main(void) 
{
  int i, id;
  double x;

#ifdef SAVE
  for (i = 0; i < 10000; i++) {
    x = rnd0();
    id = i % 10000;
    if (id <= 10)
      printf("%6d: %16.14f\n", id, x);
  }
  mtsave(NULL);
#else
  for (i = 0; i < 20000; i++) {
    x = rnd0();
    id = i % 10000;
    if (id <= 10)
      printf("%6d: %16.14f\n", id, x);
  }
#endif
  return 0;
}
