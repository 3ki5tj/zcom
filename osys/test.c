#include "osys.c"

int main(void)
{
  char *out, **fnls;
  int i, n;

  out = sysrun("df", NULL, 0);
  printf("Output:\n%s\n\n", out);
  ssdel(out);

  if ((fnls = fnglob("*.c", &n, NULL, 0)) != NULL) {
    for (i = 0; i < n; i++) {
      printf("file %d: %s\n", i+1, fnls[i]);
    }
  }
  if (fnls) {
    ssdel(fnls[0]);
    free(fnls);
  }
  return 0;
}
