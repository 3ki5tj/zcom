#include "md.h"

int main(void)
{
  rv3_t x3[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}};
  rv2_t x2[4] = {{0, 0}, {0, 1}, {.8f, 1.6f}, {1.8f, 1.6f}};
  int n = 4;

  md_shiftcom3d(x3, n);
  md_shiftcom2d(x2, n);
  return 0;
}
