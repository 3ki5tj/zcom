#include <stdio.h>
#include "rv3.h"
#include "pdb.c"

int main(int argc, char **argv)
{
  pdbmodel_t *m;
  pdbaac_t *c;
  const char *fn;

  fn = (argc > 1) ? argv[1] : "test.pdb";
  if ((m = pdbm_read(fn, 2)) == NULL) {
    return -1;
  }
  pdbm_write(m, "out.pdb");
  c = pdbaac_parse(m, 1);
  die_if (c == NULL, "bad pdb\n");
  pdbaac_free(c);
  pdbm_free(m);
  return 0;
}

