#include "blkmem.h"

typedef struct {
  int ival;
  double dval[4];
} mystruct_t;


int main(void)
{
  blkmem_t *bm;
  size_t blksz = 32, i;
  mystruct_t *p;

  bm = blkmem_open(sizeof(mystruct_t), blksz);
  die_if (bm == NULL, "no bm %u\n", (unsigned) blksz);
  for (i = 0; i < 100000u; i++) {
    p = blkmem_new(bm, i % 4 + 3);
    if (p != NULL) {
      p->ival = i + 1;
      p->dval[0] = i * i;
    } else {
      fprintf(stderr, "no more memories, i %u\n", (unsigned) i);
      break;
    }
  }
  blkmem_print(bm, "mystruct");
  blkmem_close(bm);
  return 0;
}
