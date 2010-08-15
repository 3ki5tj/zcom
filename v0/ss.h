#ifndef SS_H__
#define SS_H__
#include <string.h>

#define sscopy(s, t) sscopy_x(&(s), t)
#define sscat(s, t)  sscat_x(&(s), t)
#define ssfinish() ssmanage(0)
#define ssshrink() ssmanage(1)

char *ssnew(const char *);
void ssdelete(char *);
void ssmanage(int);
char *sscopy_x(char **, const char *);
char *sscat_x(char **, const char *);

#endif

