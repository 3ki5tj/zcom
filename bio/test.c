#include <stdio.h>
#include <math.h>
#include "bio.h"

#define datafile "io.bin"
#define N 4
int ival;
double dval;
int iarr[N];
double darr[N];

static void print(void)
{
  int i;

  printf("ivar = 0x%X, dval = %g\n", ival, dval);
  for (i = 0; i < N; i++)
    printf("%d: %10d %18.14f\n", i, iarr[i], darr[i]);
  printf("\n");
}

static int bread(const char *fname)
{
  FILE *fp;
  int err, endn, size;

  if ((fp=fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot read file %s\n", fname);
    return 1;
  }

  BIO_INIT_ENDIAN(size, sizeof(int));
  BIO_RMI(size, sizeof(double));

  BIO_RI(ival);
  BIO_RD(dval);
  BIO_RIARR(iarr, N);
  BIO_RDARR(darr, N);
  
  printf("successfully read file %s\n", fname);
  print();
  getchar();

  fclose(fp);
  return 0;
ERR:
  fclose(fp);
  return 1;
}

static int bwrite(const char *fname)
{
  FILE *fp;
  int size;

  if ((fp=fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write file %s\n", fname);
    return 1;
  }
  
  size = sizeof(int);
  BIO_WI(size);
  size = sizeof(double);
  BIO_WI(size);
  BIO_WI(ival);
  BIO_WD(dval);
  BIO_WIARR(iarr, N);
  BIO_WDARR(darr, N);

  fclose(fp);
  
  printf("successfully wrote file %s\n", fname);
  getchar();
  return 0;
ERR:
  fclose(fp);
  return 1;
}

static void initdata(void)
{
  int i, pwr;
  double x;

  ival = 0x1234ABCD;
  dval = M_PI;

  for (pwr = 1, i = 0; i < N; i++, pwr *= 2) {
    iarr[i] = pwr;
  }

  for (x = 0.6, i = 0; i < N; i++) {
    x = 3.8 * x * (1 - x);
    darr[i] = x;
  }
}

int main(void)
{
  if (bread(datafile) != 0) {
    fprintf(stderr, "no initial file\n");
  }
  initdata();
  bwrite(datafile);
  bread(datafile);
  return 0;
}

