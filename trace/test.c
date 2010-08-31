#include <stdio.h>
#include "trace.h"

int main(void)
{
  char msg[] = "Hello World, hope this is long enough! ";
  int i;


#define test(func, n, m) \
  printf("testing %s for %d repetitions ... \n", #func, n); \
  getchar(); \
  for (i = 1; i <= n; i++) { \
    func("%6d: %24.12e %s\n", i, 3.14, msg); \
    if (i % m == 0) { \
      printf("press Enter to continue, i = %d", i); \
      getchar(); \
    } \
  } \
  func(NULL)

  wtrace_buf("%@freq=",  500);
  test(wtrace,     1321, 100);
  test(wtrace_buf, 1321, 100);

  return 0;
}
