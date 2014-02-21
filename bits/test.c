#include <stdio.h>
#include <stdlib.h>

#define WORDBITS_ 32
#include "bits.h"

static void testfirst(void)
{
  int i, k;
  unsigned int x = 32143;

  for ( i = 0; i < 1000000; i++ ) {
    x = x * 22695477u + 1u;
    k = bitfirst(x);
    if ( !(x & MKBIT(k)) || (x & MKBITSMASK(k)) ) {
      fprintf(stderr, "failed test x %#x, k %d\n", x, k);
      exit(1);
    }
  }
  fprintf(stderr, "passed the bitfirst test\n");
}


static void testcount(void)
{
  uint64_t x = (uint64_t) -1;

  printf(" 8: %d\n", BITCOUNT8(x));
  printf("16: %d\n", BITCOUNT16(x));
  printf("24: %d\n", BITCOUNT24(x));
  printf("32: %d\n", BITCOUNT32(x));
  printf("40: %d\n", BITCOUNT40(x));
  printf("48: %d\n", BITCOUNT48(x));
  printf("56: %d\n", BITCOUNT56(x));
  printf("64: %d\n", BITCOUNT64(x));
}



int main(void)
{
  testfirst();
  testcount();
  return 0;
}
