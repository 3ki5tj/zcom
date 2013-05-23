#include <stdio.h>
#include <string.h>
/*
#define HAVEVAM 1
*/
#include "util.h"

int main(void)
{
  char s[] = "hello world";

  die_if (strlen(s) > 8,
      "string \"%s\" is too long\n", s);
  return 0;
}
