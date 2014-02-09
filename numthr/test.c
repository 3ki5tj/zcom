#include "numthr.c"

void test_modpsqrt(int p)
{
  int n, r;

  for (n = 1; n < p; n++) {
    r = modpsqrt(n, p);
    if ( r < 0 ) continue;
    printf("modpsqrt(%d, %d) = %d, sqr = %d\n", n, p, r, modpow(r, 2, p));
  }
}


int main(void)
{
  test_modpsqrt(43);
  return 0;
}
