#ifndef SS_H__
#define SS_H__
#include <string.h>

enum { SSCAT=1, SSDELETE=2, SSSHRINK=3, SSSINGLE=0x1000 };

#define ssnew(t)       sscpycatx(NULL, (t), 0)
#define sscpy(s, t)    sscpycatx(&(s), (t), 0)
#define sscat(s, t)    sscpycatx(&(s), (t), SSCAT)
#define ssdel(s)       ssmanage((s), SSDELETE|SSSINGLE)
#define ssshr(s)       ssmanage((s), SSSHRINK|SSSINGLE)
#define ssdelall()     ssmanage(NULL, SSDELETE)
#define ssshrall()     ssmanage(NULL, SSHRINK)
#define ssfgets(s, pn, fp)    ssfgetx(&(s), (pn), '\n', (fp))
#define ssfgetall(s, pn, fp)  ssfgetx(&(s), (pn), EOF, (fp))

void ssmanage(char *, unsigned);
char *sscpycatx(char **, const char *, unsigned);
char *ssfgetx(char **, size_t *, int, FILE *fp);

#endif

