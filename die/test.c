#include <stdio.h>
#include <string.h>
#include "die.h"

int main(void)
{
  char s[] = "hello world";

  die_if (strlen(s) > 8, 
      "string is too long\n");
  return 0;
}
