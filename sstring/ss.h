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
#define ssfgets(s, pn, fp)    ssfgetx(&(s), pn, '\n', fp)
#define ssfgetall(s, pn, fp)  ssfgetx(&(s), pn, EOF, fp)

void ssmanage(char *, unsigned);
char *sscpyx(char **, const char *);
char *sscatx(char **, const char *);
char *ssfgetx(char **, size_t *, int, FILE *fp);

#endif

