#include <stdio.h>
#include <string.h>
#include "util.h"

int main(void)
{
  char fn[] = "util.h";
  char fa[] = "../zcom.h", fb[] = "zcom.bak";

  printf("%s exists? %s\n", fn, (fexists(fn) ? "yes" : "no"));
  if ( fexists(fa) ) {
    copyfile(fa, fb);
    printf("copy %s to %s\n", fa, fb);
  }
  return 0;
}
