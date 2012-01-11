/* test loading previous configuration */
#include "clus.c"

int main(void)
{
  clsys_t *cls;
  int verbose = 2;

  cls = cls_read("clus.txt", NULL, NULL);
  cls = cls_zalgo(cls, 10*cls->np*cls->np, CLUS_HEATBATH, 0.2, 20.0, 100, 1000, verbose);
  cls_write(cls, "clus1.txt", NULL, NULL, 0);
  cls_free(cls, 1);

  return 0;
}
