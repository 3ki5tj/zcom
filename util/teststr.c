#include "util.h"



static void test_substr(void)
{
  const char *s = "Hello, world";
  char t[5];
  int start = 2, len = 4;

  printf("s = [%s]\n", s);
/*
  substr() is currently removed for potential name conflict
  substr(t, s, start, len);
*/
  strcpy_sf(t, s + start, len + 1);
  printf("substr(t, s, %d, %d) --> t = [%s]\n",
      start, len, t);
}



static void test_strcmpnc(void)
{
  const char *s = "Hello", *t = "HELLo";
  const char *s2 = "Abcd", *t2 = "aBCdF";

  printf("\nTesting strcmpnc...\n");
  printf("s = [%s], t = [%s]\n", s, t);
  printf("strcmp(s,t) == %d, strcmpnc(s,t) == %d\n",
      strcmp(s, t), strcmpnc(s, t));
  printf("s2 = [%s], t2 = [%s]\n", s2, t2);
  printf("strncmp(s2,t2,4) == %d, strncmpnc(s2,t2,4) == %d\n",
      strncmp(s2, t2, 4), strncmpnc(s2, t2, 4));
  printf("\n");
}



static void test_strip(void)
{
  char s[128];
  const char *src = "  Hello world!\t ";

  printf("\nTesting strip(s)...\n");
  printf("s = [%s]\n", src);
  printf("lstrip(s) = [%s]\n", lstrip(strcpy(s, src)));
  printf("rstrip(s) = [%s]\n", rstrip(strcpy(s, src)));
  printf("strip(s) = [%s]\n", strip(strcpy(s, src)));
  printf("\n");
}



int main(void)
{
#ifdef __unix__
  printf("unix operating system\n");
#endif
#ifdef __WIN32__
  printf("unix operating system\n");
#endif
  test_strcmpnc();
  test_strip();
  test_substr();
  return 0;
}

