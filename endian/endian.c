#ifndef ENDIAN_C__
#define ENDIAN_C__

#include "endian.h"

/* Correct endianness from the current operating system to the desired one
 * operates on an array of `n' objects, each of `size' bytes
 * The target endianness is specified by endtar, whereas
 * the endianness of the current operating system is automatically computed.
 * A conversion takes place only if the two differ.
 * If output is not NULL, the enidan-corrected array is saved there, 
 * otherwise an in-place conversion is assumed and input is returned
 * */
void *fix_endian(void *output, void *input, size_t size, size_t n, int endtar)
{
  static int endsys = -1;
  unsigned char *dest = (unsigned char *) output;
  unsigned char *src  = (unsigned char *) input;
  unsigned char *p, *q, ch;
  size_t i, r, j;

  if (endsys < 0) { /* initial determine the machine's endianess */
    unsigned feff = 0xFEFF; /* assume unsigned is at least 16-bit */
      
    p  = (unsigned char *) &feff;
    endsys = (*p == 0xFF) ? 0 : 1;
#ifdef ENDIAN_DBG_
    fprintf(stderr, "The current system is %s-endian.\n", 
        endsys ? "big" : "little");
#endif
  }

  if (endtar == endsys) { /* no conversion is required */
    /* out-of-place: just copy array; in-place: do nothing */
    return (dest != NULL) ? memcpy(dest, src, size * n) : src;
  } 

  for (p = src, q = dest, j = 0; j < n; j++, p += size) {
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

#endif

