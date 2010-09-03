#include <stdio.h>
#include "endn.h"

static void pint(int i)
{
  int n;
  unsigned char *p;

  printf("0x%08X,%12d: ", i, i);
  p = (unsigned char*)&i;
  for (n = 0; n < sizeof(int);  n++) {
    printf( "0x%02X ", p[n]);
  }
  printf("--> \n");
}

/* in-place conversion */
#define END2BIG(from) endn_converti(&(from), sizeof(from), 1, 1)
#define END2LIT(from) endn_converti(&(from), sizeof(from), 1, 0)

/* out-of-place conversion */
#define END2BIGO(to, from) endn_converti(&(to), &(from), sizeof(from), 1, 1)
#define END2LITO(to, from) endn_converti(&(to), &(from), sizeof(from), 1, 0)

static void test1(void)
{
  int i = 0xFEFF, j;

  printf("\n\nTesting in-place ...\n");
  i = 0xFEFF;
  printf("i before END2BIG\n");
  pint(i);
  END2BIG(i);
  printf("i after  END2BIG\n");
  pint(i);
  END2BIG(i);
  printf("i after  END2BIG again\n");
  pint(i);

  i = 0xFEFF;
  printf("\n\ni before END2LIT\n");
  pint(i);
  END2LIT(i);
  printf("i after  END2LIT\n");
  pint(i);

  printf("\n\nTesting out-of-place ...\n");
  i = 0xFEFF;
  j = 0;
  printf("\n\ni, j before END2BIGO(j, i)\n");
  pint(i);
  pint(j);
  END2BIGO(j, i);
  printf("i, j after  END2BIGO(j, i)\n");
  pint(i);
  pint(j);

  i = 0xFEFF;
  j = 0;
  printf("\n\ni, j before END2LITO(j, i)\n");
  pint(i);
  pint(j);
  END2LITO(j, i);
  printf("i, j after  END2LITO(j, i)\n");
  pint(i);
  pint(j);
  
  printf("press any key to continue..."); getchar();
}

#define N 3
static void parr(unsigned s[], const char *msg)
{
  int i;
  unsigned char *p;
  
  printf("%s:\n", msg);
  for (i = 0; i < N; i++)
    printf("0x%08X ", s[i]);

  printf(" | ");
  p = (unsigned char*)s;
  for (i = 0; i < N*sizeof(int);  i++) 
    printf( "0x%02X ", p[i]);
  printf("\n");
  fflush(stdout);
}

/* test array case */
static void test2(void)
{
  unsigned src[N] = {0xFEFF, 0x123456, 0x7890ABCD};
  unsigned s[N], t[N];

  printf("\n\nTesting in-place ...\n\n");
  memcpy(s, src, N*sizeof(s[0]));
  parr(s, "before to big endian");
  fix_endian_inp(s, sizeof(s[0]), N, 1); /* big endian */
  parr(s, "after to big endian");

  memcpy(s, src, N*sizeof(s[0]));
  parr(s, "before to little endian");
  fix_endian_inp(s, sizeof(s[0]), N, 0);
  parr(s, "after to little endian");

  printf("\n\nTesting out-of-place ...\n\n");
  memcpy(s, src,  N*sizeof(s[0]));
  memset(t, '\0', N*sizeof(s[0]));
  parr(s, "s: before to big endian");
  parr(t, "t: before to big endian");
  fix_endian(t, s, sizeof(s[0]), N, 1); /* big endian */
  parr(s, "s: after to big endian");
  parr(t, "t: after to big endian");

  memcpy(s, src,  N*sizeof(s[0]));
  memset(t, '\0', N*sizeof(s[0]));
  parr(s, "s: before to little endian");
  parr(t, "t: before to little endian");
  fix_endian(t, s, sizeof(s[0]), N, 0);
  parr(s, "s: after to little endian");
  parr(t, "t: after to little endian");

  printf("press any key to continue..."); getchar();
}

int main(void)
{
  test1();
  test2();
  return 0;
}

