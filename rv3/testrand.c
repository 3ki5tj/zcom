#include "rv3.h"
#include <time.h>



static void speed_dir(int n)
{
  int i;
  clock_t t;
  real v[3];

  t = clock();
  for (i = 0; i < n; i++)
    rv3_randdir0(v);
  printf("randdir0() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++) {
    while ( rv3_sqr(rv3_randunif(v, -1, 1)) >= 1 ) ;
    rv3_normalize(v);
  }
  printf("randdir0b() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++) {
    rv3_randgaus(v);
    rv3_normalize(v);
  }
  printf("randdir0c() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);
}



static void speed_ball(int n)
{
  int i;
  clock_t t;
  real v[3];

  t = clock();
  for (i = 0; i < n; i++)
    rv3_randball0(v);
  printf("randball0() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++) {
    real r = (real) pow(rand01(), 1./3);
    rv3_randdir(v, r);
  }
  printf("randball0b() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);
}



static void print(int n)
{
  FILE *fp;
  int i;
  rv3_t v;

  fp = fopen("a.dat", "w");
  for (i = 0; i < n; i++) {
    rv3_randball0(v);
    fprintf(fp, "%g %g %g\n", v[0], v[1], v[2]);
  }
  fclose(fp);

  fp = fopen("b.dat", "w");
  for (i = 0; i < n; i++) {
    real r = (real) pow(rand01(), 1./3);
    rv3_randdir(v, r);
    fprintf(fp, "%g %g %g\n", v[0], v[1], v[2]);
  }
  fclose(fp);
}



int main(void)
{
  speed_dir(20000000);
  speed_ball(20000000);
  print(10000);
  return 0;
}
