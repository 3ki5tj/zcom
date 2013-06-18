#define HAVEVAM
#define NEED_WTRACE
#include "log.h"

int main(void)
{
  char msg[] = "Hello World, hope this is long enough! ";
  int i;


#define test(func, n) \
  printf("testing %s for %d repetitions ... \n", #func, n); \
  for (i = 1; i <= n; i++) { \
    func("%6d: %24.12e %s\n", i, 3.14, msg); \
  } \
  func(NULL, NULL)

  wtrace("%@freq=",  2);
  test(wtrace, 1);

  return 0;
}
