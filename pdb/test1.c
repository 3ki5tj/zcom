#include <stdio.h>
#include "rv3.h"
#include "pdb.h"

int main(void)
{
  pdbmodel_t *m;
  pdbaabb_t *bb;

  if ((m = pdbload0("test.pdb", 2)) == NULL) {
    return -1;
  }
  bb = pdbgetaabb(m, 1);
  return 0;
}

