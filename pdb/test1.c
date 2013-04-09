#include <stdio.h>
#include "rv3.h"
#include "pdb.c"

int main(int argc, char **argv)
{
  pdbmodel_t *m;
  pdbaac_t *c;
  const char *fn;
  int *se, i, ng;

  fn = (argc > 1) ? argv[1] : "test.pdb";
  die_if ((m = pdbm_read(fn, 2)) == NULL, "cannot read %s\n", fn);
  pdbm_write(m, "out.pdb");
  c = pdbaac_parse(m, 1);
  die_if (c == NULL, "bad pdb\n");
  
  ng = pdbaac_parsehelices(c, &se);
  for (i = 0; i < ng; i++) /* print helix groups */
    printf("Helix %2d: %d - %d\n", i + 1, se[2*i] + 1, se[2*i+1]);

  free(se);
  pdbaac_free(c);
  pdbm_free(m);
  return 0;
}

