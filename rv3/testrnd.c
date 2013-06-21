#include "rv3.h"
#include <time.h>


/* compare the speed of rndball0(v) and the alternative implementaion:
 * rnddir(v, rnd0()). Although the later is slightly faster, the former is simpler */
static void speed(int n)
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
    rv3_rnddir(v, rnd0());
  printf("rnddir() used %gs\n", (double)(clock() - t)/CLOCKS_PER_SEC);
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
    rv3_rnddir(v, rnd0());
    fprintf(fp, "%g %g %g\n", v[0], v[1], v[2]);
  }
  fclose(fp);
}



int main(void)
{
  speed(20000000);
  print(10000);
  return 0;
}
