#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef SS_H__
#define SS_H__
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

/* operation code */
enum { SSCAT = 1, SSDELETE = 2, SSSHRINK = 3, SSSINGLE = 0x1000 };

#define ssnew(n)       sscpycatx(NULL, NULL, (n),    0)
#define ssdup(t)       sscpycatx(NULL, (t),   0,     0)
#define sscpy(s, t)    sscpycatx(&(s), (t),   0,     0)
#define sscat(s, t)    sscpycatx(&(s), (t),   0, SSCAT)
#define ssdel(s)       ssmanage((s), SSDELETE|SSSINGLE)
#define ssdelete(s)    { ssdel(s); (s)=NULL; }
#define ssshrink(s)    ssmanage((s), SSSHRINK|SSSINGLE)
#define ssdelall()     ssmanage(NULL, SSDELETE)
#define ssshrinkall()  ssmanage(NULL, SSHRINK)
#define ssfgets(s, pn, fp)    ssfgetx(&(s), (pn), '\n', (fp))
#define ssfgetall(s, pn, fp)  ssfgetx(&(s), (pn), EOF, (fp))

INLINE int   ssmanage(char *, unsigned);
INLINE char *sscpycatx(char **, const char *, size_t, unsigned);
INLINE char *ssfgetx(char **, size_t *, int, FILE *fp);

INLINE char **ssparse(char *s, int *pn, const char *delim);

#endif

