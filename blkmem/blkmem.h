#include "util.h"
#ifndef BLKMEM_H__
#define BLKMEM_H__



/* block memory allocator */

typedef struct tagblkmemlls_t {
  char *ptr; /* allocated memory */
  size_t n; /* number of item used */
  struct tagblkmemlls_t *next;
} blkmemlls_t;

typedef struct {
  size_t objsz; /* memory for a single object */
  size_t blksz; /* number of objects to allocate each time */
  blkmemlls_t *ls, *last;
  int err; /* we are out of memory */
} blkmem_t;



INLINE blkmem_t *blkmem_open(size_t objsz, size_t blksz)
{
  blkmem_t *bm;

  if ((bm = calloc(1, sizeof(*bm))) == NULL) {
    fprintf(stderr, "blkmem: no memory for blkmem_t\n");
    return NULL;
  }
  bm->objsz = objsz;
  bm->blksz = blksz;
  bm->ls = NULL;
  bm->last = NULL; /* for quick look up */
  bm->err = 0;
  return bm;
}



INLINE void blkmem_close(blkmem_t *bm)
{
  blkmemlls_t *ls, *ls1;

  for (ls = bm->ls; ls; ls = ls1) {
    ls1 = ls->next;
    free(ls->ptr);
    free(ls);
  }
  free(bm);
}



/* allocate an empty list */
INLINE blkmemlls_t *blkmem_newlls(size_t objsz, size_t blksz, int *err)
{
  blkmemlls_t *ls;

  if (*err) return NULL;
  if ((ls = calloc(1, sizeof(*ls))) == NULL) {
    fprintf(stderr, "blkmem: no memory for the first list\n");
    *err = 1;
    return NULL;
  }
  ls->n = 0;
  ls->next = NULL;
  if ((ls->ptr = calloc(blksz, objsz)) == NULL) {
    fprintf(stderr, "blkmem: no memory for list objects %u x %u\n",
        (unsigned) blksz, (unsigned) objsz);
    *err = 1;
    return NULL;
  }
  return ls;
}




#define blkmem_new1(bm) blkmem_newlow(bm, 1, 1)
#define blkmem_new(bm, k) blkmem_newlow(bm, k, 1)

/* allocate k new objects of size bm->objsz
 * if `fromlast' is true, we do not search for holes in previous lists
 * otherwise we search from the beginning of the list, which is much slower */
INLINE void *blkmem_newlow(blkmem_t *bm, size_t k, int fromlast)
{
  blkmemlls_t *ls;
  size_t blksz = bm->blksz, objsz = bm->objsz, n;

  if (k == 0) {
    return NULL;
  } else if (k >= blksz) {
    fprintf(stderr, "blkmem: size %u > %u\n", (unsigned) k, (unsigned) blksz);
    return NULL;
  }

  if (bm->ls == NULL) { /* initial empty list */
    if ((ls = blkmem_newlls(objsz, blksz, &bm->err)) == NULL)
      return NULL;
    bm->ls = bm->last = ls;
    //printf("new list\n"); getchar();
  } else {
    if (fromlast) /* search from the last item */
      ls = bm->last;
    else /* search from the beginning */
      ls = bm->ls;
    //die_if (ls == NULL, "blkmem, last list is zero k %d\n", k);
    /* loop to the end of the list
     * or stop if there are enough spaces left */
    while (1) {
      if (ls->n + k <= blksz) { /* find a hole */
        //printf("current list, ls->n %d\n", ls->n); getchar();
        break;
      } else if (ls->next == NULL) { /* end of the list, extend it */
        if ((ls->next = blkmem_newlls(blksz, objsz, &bm->err)) == NULL)
          return NULL;
        ls = ls->next;
        bm->last = ls;
        //printf("extended list, ls->n %d\n", ls->n); getchar();
        break;
      } else { /* check the next list */
        ls = ls->next;
      }
    }
  }
  /* we still have memory */
  n = ls->n;
  ls->n = n + k;
  return (void *) (ls->ptr + n * objsz);
}



#define blkmem_print(bm, name) blkmem_fprint(bm, name, stdout)

/* print usage */
INLINE void blkmem_fprint(const blkmem_t *bm, const char *name, FILE *fp)
{
  blkmemlls_t *ls = bm->ls;

  fprintf(fp, "(blkmem)%s: %u (objsz) x %u (blksz): ",
      name ? name : "", (unsigned) bm->objsz, (unsigned) bm->blksz);
  for (; ls; ls = ls->next)
    fprintf(fp, "%u ", (unsigned) ls->n);
  fprintf(fp, "\n");
}


#endif /* BLKMEM_H__ */
