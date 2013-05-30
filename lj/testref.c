#include "lj.c"

int main(void)
{
  double G, F, U, P, rho, T = 1.6;
  for (rho = 0.1; rho <= 1.01; rho += 0.1) {
    U = lj_eos3d(rho, T, &P, &F, &G);
    printf("rho %g, T %g, P %g, U %g, F %g, G %g\n", rho, T, P, U, F, G);
  }
  return 0;
}
