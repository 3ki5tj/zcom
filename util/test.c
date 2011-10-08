#include <stdio.h>
#include <string.h>
/*
#define HAVEVAM 1
*/
#include "util.c"

int main(void)
{
  char s[] = "hello world";

#ifdef USE_MSG_IF
  msg_if (strcmp(s, "Hello World") != 0,
      "strings are different\n");
#endif
  die_if (strlen(s) > 8, 
      "string is too long\n");
  return 0;
}
