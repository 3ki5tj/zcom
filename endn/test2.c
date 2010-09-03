/* convert an integer array */
#include <stdio.h>
#include "endn.h"

#define N 3
unsigned src[N] = {0xFEFF, 0x123456, 0x7890ABCD};
unsigned s[N], t[N];
char *endname[2] = {"little endian", "big endian"};

static void parr(unsigned s[], const char *msg, int endn)
{
  int i;
  unsigned char *p;
  
  printf("%s %s:\n", msg, endname[endn]);
  for (i = 0; i < N; i++)
    printf("0x%08X ", s[i]);

  printf(" | ");
  p = (unsigned char*)s;
  for (i = 0; i < N*sizeof(int);  i++) 
    printf( "0x%02X ", p[i]);
  printf("\n");
  fflush(stdout);
}

static void inplace(int endn)
{
  memcpy(s, src, N*sizeof(s[0]));
  parr(s, "before to", endn);
  endn_converti(s, sizeof(s[0]), N, endn);
  parr(s, "after to", endn);
}

static void outofplace(int endn)
{
  memcpy(s, src,  N*sizeof(s[0]));
  memset(t, '\0', N*sizeof(s[0]));
  parr(s, "s: before to", endn);
  parr(t, "t: before to", endn);
  endn_convert(t, s, sizeof(s[0]), N, endn); /* big endian */
  parr(s, "s: after to", endn);
  parr(t, "t: after to", endn);
}

/* test array case */
int main(void)
{
  printf("\n\nTesting in-place ...\n\n");
  inplace(1); /* big endian */
  inplace(0); /* little endian */

  printf("\n\nTesting out-of-place ...\n\n");
  outofplace(1); /* big endian */
  outofplace(0); /* little endian */

  return 0;
}

