#ifndef SS_H__
#define SS_H__
#include <string.h>

enum{ SSDELETE, SSSHRINK, SSSINGLE=0x1000};

#define ssnew(t)       sscpyx(NULL, t)
#define sscpy(s, t)    sscpyx(&(s), t)
#define sscat(s, t)    sscatx(&(s), t)
#define ssdel(s)       ssmanage(s, SSDELETE|SSSINGLE)
#define ssshr(s)       ssmanage(s, SSSHRINK|SSSINGLE)
#define ssdelall()     ssmanage(NULL, SSDELETE)
#define ssshrall()     ssmanage(NULL, SSHRINK)

void ssmanage(char *, unsigned);
char *sscpyx(char **, const char *);
char *sscatx(char **, const char *);

#endif

