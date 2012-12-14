/* print transition matrices according to code 
 * Copyright (C) 2012 Cheng Zhang */
#include "spgr.h"

int main(int argc, char **argv)
{
  int rot[SPGR_OPMAX][3][3];
  double trans[SPGR_OPMAX][3];
  int i, cnt, hall;
  char *code = "P4_122";

  if (argc > 1) code = argv[1];

  cnt = spgr_getmatstr(code, &hall, rot, trans);
  printf("# %d %d %d %s\n", cnt,
      spacegroup_types[hall].number, hall,
      spgr_getpgname(hall, 1));
  for (i = 0; i < cnt; i++) {
    printf("%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f "
        "%8.5f %8.5f %8.5f\n",
        1.*rot[i][0][0], 1.*rot[i][0][1], 1.*rot[i][0][2],
        1.*rot[i][1][0], 1.*rot[i][1][1], 1.*rot[i][1][2],
        1.*rot[i][2][0], 1.*rot[i][2][1], 1.*rot[i][2][2],
        trans[i][0], trans[i][1], trans[i][2]);
  }
  return 0;
}
