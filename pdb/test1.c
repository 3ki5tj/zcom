#include <stdio.h>
#include "rv3.h"
#include "pdb.c"

int main(void)
{
  pdbmodel_t *m;
  pdbaac_t *c;

  if ((m = pdbm_read("test.pdb", 2)) == NULL) {
    return -1;
  }
  pdbm_write(m, "out.pdb");
  c = pdbaac_parse(m, 1);
  pdbaac_free(c);
  pdbm_free(m);
  return 0;
}

