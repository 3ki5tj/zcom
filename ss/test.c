#include <stdio.h>
#include "ss.h"

static int test_ssfgets(void)
{
  char *s = NULL;
  const char *fname = "ss.c";
  FILE *fp;
  int i;
  size_t n = 0;

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot open file, %s\n", fname);
    return -1;
  }
  for (i = 1; ssfgets(s, &n, fp); i++)
    printf("%4d %4u | %s", i, n, s);
  rewind(fp);
  ssfgetall(s, &n, fp);
  printf("File %s has %u characters\n", fname, n);
  fclose(fp);
  ssdel(s);
  return 0;
}



int main(void)
{
  int i;
  const char *sa="A", *sb="B";
  char *s = NULL, *t = NULL, *r = NULL;

  r = ssnew(512);
  for (i = 0; i < 512; i++)
    r[i] = 'j';
  r[i] = '\0';

  s = ssdup(sa);
  t = ssdup(sb);
  for (i = 1; i <= 16; i++) {
    sscpy(r, s); /* make a copy of the old s */
    sscat(s, t); /* s = s+t */
    sscpy(t, r); /* t = old s */
    printf("%3d: %6u, %s\n", i, strlen(s), s);
  }
  sscpy(s, "new");
  ssshrink(s);
  ssdelall();

  test_ssfgets();
  return 0;
}

