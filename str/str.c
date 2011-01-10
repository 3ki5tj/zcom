#ifndef STR_C__
#define STR_C__

#include "str.h"

/* convert strings to different cases.
 * 0 for success,
 * 1 for lack of space, or NULL pointers.
 * Unlike the case of strncpy(), s is always null-terminated.
 * */
char *strcnv(char *s, const char *t, size_t size_s, unsigned flags)
{
  size_t i, j;
  int cs, ct, docase, doupper;

  docase  = flags & ZSTR_CASE; /* do case conversion */
  doupper = flags & ZSTR_UPPER;

  if (size_s == 0 || s == NULL || t == NULL) return s;
  /* t[size_s-1] should be '\0' */
  i = 0;
  if (flags & ZSTR_CAT) while(s[i]) i++;
  for (j = 0; i < size_s-1; i++, j++) {
    ct = t[j];
    if (docase && (ct != 0)) {
      ct = (unsigned char)(char)ct;
      if (doupper) {
        cs = toupper(ct);
      } else {
        cs = tolower(ct);
      }
    } else {
      cs = ct;
    }
    s[i] = (char)cs;
    if (ct == 0) break;
  }
  if (i == size_s-1) s[i] = '\0';
  return s;
}

/* compare strings ignoring cases */
int strncmpnc(const char *s, const char *t, int n)
{
  int i, cs, ct;

  if (s == NULL || t == NULL) return 0;
  for (i = 0; ; i++) {
    if (i >= n) return 0;
    cs = s[i];
    ct = t[i];
    if (cs == 0 || ct == 0) break;
    cs = tolower( (unsigned char) (int) cs );
    ct = tolower( (unsigned char) (int) ct );
    if (cs != ct) break;
  }
  return cs-ct;
}

/* remove leading and trailing spaces */
char *stripx(char *s, unsigned flags)
{
  char *p;
  /* strip the leading */
  if (flags & ZSTR_LEFT) {
    for (p = s; isspace(*p); p++) ;
    if (*p == '\0') {
      *s = '\0';
    } else memmove(s, p, strlen(p)+1);
  }

  /* strip the ending */
  if (flags & ZSTR_RIGHT) {
    for (p = s + strlen(s) - 1; p >= s && isspace(*p); p--)
      *p = '\0';
  }
  return s;
}

#endif
