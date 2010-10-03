#ifndef ENDN_C__
#define ENDN_C__

#include "endn.h"

/* return the system endian, 1: big endian, 0: little endian */
int endn_system(void)
{
  unsigned feff = 0xFEFF; /* assume unsigned is at least 16-bit */
  unsigned char *p;
      
  p  = (unsigned char *) &feff;
  return (*p == 0xFF) ? 0 : 1;
}

/* change endianness for items of `size' in src 
 * results are saved to dest, if it's not NULL, 
 * or an in-place conversion is performed  */
void *endn_flip(void *dest, void *src, size_t size, size_t n)
{
  unsigned char *p, *q, ch;
  size_t i, r, j;

  p = (unsigned char *) src;
  q = (unsigned char *) dest;
  for (j = 0; j < n; j++, p += size) {
    /* reverse bytes for each object */
    if (dest == NULL) { /* in-place conversion */
      for (i = 0; i < (size/2); i++) {
        r = size - i - 1;
        ch   = p[i];
        p[i] = p[r];
        p[r] = ch;
      }
    } else { /* out-of-place conversion */
      for (i = 0; i < size; i++)
        q[size - i - 1] = p[i];
      q += size;
    }
  }
  return (dest != NULL) ? dest : src;  
}

#ifdef ENDN_LEGACY
/* correct endianness from the current operating system to the desired one
 * operates on an array of `n' objects, each of `size' bytes
 * The target endianness is specified by `endn', whereas
 * the endianness of the current operating system is automatically computed.
 * A conversion takes place only if the two differ.
 * If dest is not NULL, the enidan-corrected array is saved there, 
 * otherwise an in-place conversion is assumed and src is returned
 * */
void *endn_convert(void *dest, void *src, size_t size, size_t n, int endn)
{
  static int endsys = -1;

  /* initial determine the machine's endianess */
  if (endsys < 0) endsys = endn_system();

  if (endn == endsys) { /* no conversion is required */
    /* out-of-place: just copy array; in-place: do nothing */
    return (dest != NULL) ? memcpy(dest, src, size * n) : src;
  } else {
    return endn_flip(dest, src, size, n);
  }
}
#endif

/* write data in ptr to file with a specific endian 'endn' */
size_t endn_fwrite(const void *ptr, size_t size, size_t n, FILE *fp, int endn)
{
  static int endsys = -1;
  size_t i;
  unsigned char buf[16], *p;

  /* initial determine the machine's endianess */
  if (endsys < 0) endsys = endn_system();
  if (endn == endsys)
    return fwrite(ptr, size, n, fp);
  
  if (size > sizeof(buf)) {
    fprintf(stderr, "cannot perform conversion, object too large: "
        "size = %u\n", (unsigned) size);
    return 0;
  }

  p = (unsigned char *) ptr;
  for (i = 0; i < n; i++, p += size) {
    endn_flip(buf, (void *) p, size, 1);
    if (1 != fwrite(buf, size, 1, fp)) {
      fprintf(stderr, "error occurs when writing %u / %u object\n",
          (unsigned) i, (unsigned) n);
      return i;
    }
  }
  return n;
}

/* read an object test object *src, compared with *ref
 * return 0 if they are identical without endian change 
 * return 1 if changing the endianness of *src matches *ref
 * otherwise return -1 */
int endn_rmatch(void *src, const void *ref, size_t size, FILE *fp)
{
  if (1 != fread(src, size, 1, fp))
    return -1;
#ifdef ENDN_DBG
  if (size == sizeof(int))
    printf("A: 0x%X vs. 0x%X size = %u, cmp = %d\n",
      *(int *)src, *(int *)ref, (unsigned)size, 
      memcmp(src, ref, size));
#endif
  if (memcmp(src, ref,  size) == 0)
    return 0;
  /* alter the endianness, and test again */
  endn_flip(NULL, src, size, 1);
#ifdef ENDN_DBG
  if (size == sizeof(int))
    printf("B: 0x%X vs. 0x%X size = %u, cmp = %d\n", 
      *(int *)src, *(int *)ref, (unsigned)size, 
      memcmp(src, ref, size));
#endif
  return (memcmp(src, ref, size) == 0) ? 1 : -1;
}

/* special case of endn_rmatchi for integer, convenient for iref 
 * could be something like sizeof(int), which has no address */
int endn_rmatchi(int *src, int iref, FILE *fp)
{
  return endn_rmatch(src, &iref, sizeof(int), fp);
}

/* read data from file to ptr with endianness changed if 'flip' is 1 
 * flip can be initialized by calling endn_rmatch() for a test object */
size_t endn_fread(void *ptr, size_t size, size_t n, FILE *fp, int flip)
{
  n = fread(ptr, size, n, fp);
  if (flip) endn_flip(NULL, ptr, size, n);
  return n;
}

#endif

