#include "lp.c"
static void test1(void)
{
  lpcore_t *lp = lpcore_open(3, 3);
  double a[4][4] = { {0, 3, 1, 2}, /* target function */
    {30, 1, 1, 3}, 
    {24, 2, 2, 5}, 
    {36, 4, 1, 2},
  };
  int j;
  lpcore_setc(lp, a[0]);
  for (j = 1; j <= 3; j++) lpcore_seta(lp, a[j], j);
  lpcore_printf(lp, "%8.3f");
  lpcore_simplex(lp);
  lpcore_printf(lp, "%8.3f");
  lpcore_close(lp);
}

static void test2(void)
{
  lpcore_t *lp = lpcore_open(2, 2);
  double a[3][3] = { {0, 2, -1}, /* target function */
    {2, 2, -1}, 
    {-4, 1, -5}, 
  };
  int j;
  lpcore_setc(lp, a[0]);
  for (j = 1; j <= 2; j++) lpcore_seta(lp, a[j], j);
  lpcore_printf(lp, "%8.3f");
  lpcore_simplex(lp);
  lpcore_printf(lp, "%8.3f");
  lpcore_close(lp);
}

int main(void)
{
  test1();
  test2();
  return 0;
}

