#include <stdio.h>
#include "ss.h"

int main(void)
{
  int i;
  char *sa="A", *sb="B";
  char *s=NULL, *t=NULL, *r=NULL;
  
  s=ssnew(sa);
  t=ssnew(sb);
  for(i=1; i<=16; i++){
    sscpy(r, s);
    sscat(s, t);
    sscpy(t, r);
    printf("%3d: %6u, %s\n", i, strlen(s), s);
  }
  sscpy(s, "new");
  ssshrink();
  ssfinish();
  return 0;
}

