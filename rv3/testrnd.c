#include "rv3.h"
#include <time.h>



static void speed_dir(int n)
{
  int i;
  clock_t t;
  real v[3];

  t = clock();
  for (i = 0; i < n; i++)
    rv3_rnddir0(v);
  printf("rnddir0() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++) {
    while ( rv3_sqr(rv3_rnd(v, -1, 1)) >= 1 ) ;
    rv3_normalize(v);
  }
  printf("rnddir0b() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++) {
    rv3_grand0(v);
    rv3_normalize(v);
  }
  printf("rnddir0c() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);
}



static void speed_ball(int n)
{
  int i;
  clock_t t;
  real v[3];

  t = clock();
  for (i = 0; i < n; i++)
    rv3_rndball0(v);
  printf("rndball0() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);

  t = clock();
  for (i = 0; i < n; i++)
    rv3_rnddir(v, pow(rnd0(), 1./3));
  printf("rndball0b() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);
}



static void print(int n)
{
  FILE *fp;
  int i;
  rv3_t v;

  fp = fopen("a.dat", "w");
  for (i = 0; i < n; i++) {
    rv3_rndball0(v);
    fprintf(fp, "%g %g %g\n", v[0], v[1], v[2]);
  }
  fclose(fp);

  fp = fopen("b.dat", "w");
  for (i = 0; i < n; i++) {
    rv3_rnddir(v, pow(rnd0(), 1./3));
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
