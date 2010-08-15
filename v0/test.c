#include <stdio.h>
#include "ss.h"

int main(void)
{
  int i,aprev=0;
  char *sa="A", *sb="B";
  char *s=NULL, *t=NULL, *r=NULL;
  
  s=ssnew(sa);
  t=ssnew(sb);
  for(i=1; i<=20; i++){
    sscopy(r, s);
    sscat(s, t);
    sscopy(t, r);
    printf("%3d: %s\n", i, s);
  }
  sscopy(s, "new");
  ssshrink();
  ssfinish();
  return 0;
}

