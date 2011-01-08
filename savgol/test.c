#include "savgol.h"

int main(void)
{
  double *c, *c2;
  int w = 2, ord = 2;

  c = savgol(w, ord, 0, 1);
  printf("\n\n");
  c2 = savgol2d(w, w, ord, 0, 1);
  printf("\n");
  free(c);
  free(c2);
  return 0;
}
