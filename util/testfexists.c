#include <stdio.h>
#include <string.h>
#include "util.h"

int main(void)
{
  char fn[] = "test.c";

  printf("%s exists? %s\n", fn, (fexists(fn) ? "yes" : "no"));
  return 0;
}
