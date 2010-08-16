#ifndef SS_H__
#define SS_H__
#include <string.h>

#define ssnew(t)    sscpyx(NULL, t)
#define sscpy(s, t) sscpyx(&(s), t)
#define sscat(s, t) sscatx(&(s), t)
#define ssfinish()  ssmanage(0)
#define ssshrink()  ssmanage(1)

void ssdelete(char *);
void ssmanage(int);
char *sscpyx(char **, const char *);
char *sscatx(char **, const char *);

#endif

