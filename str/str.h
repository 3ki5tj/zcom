#ifndef STR_H__
#define STR_H__

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#define ZSTR_LEFT   0x0001
#define ZSTR_RIGHT  0x0002
#define ZSTR_COPY   0x0004
#define ZSTR_CAT    0x0008
#define ZSTR_CASE   0x0100
#define ZSTR_UPPER  0x0200

/* copy the string and convert it to upper/lower case */
#define strcpy2u(s, t, size_s) strcnv(s, t, size_s, ZSTR_COPY|ZSTR_CASE|ZSTR_UPPER)
#define strcpy2l(s, t, size_s) strcnv(s, t, size_s, ZSTR_COPY|ZSTR_CASE)
#define strcpy_safe(s, t, size_s)   strcnv(s, t, size_s, ZSTR_COPY)
/* concatenate strings, the last parameter is the buffer size of s,
 * unlike strncat(), in which it's the number of characters from *t* to be copied.  */
#define strcat_safe(s, t, size_s)   strcnv(s, t, size_s, ZSTR_CAT)
#define strip(s)  stripx(s, ZSTR_LEFT|ZSTR_RIGHT)
#define lstrip(s) stripx(s, ZSTR_LEFT)
#define rstrip(s) stripx(s, ZSTR_RIGHT)

char *strcnv(char *s, const char *t, size_t size_s, unsigned flags);
int strcmpnc(const char *s, const char *t);
char *stripx(char *s, unsigned flags);

#endif

