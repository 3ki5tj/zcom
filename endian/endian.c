#ifndef ENDIAN_C__
#define ENDIAN_C__

#include "endian.h"

/* Correct endianness from the current operating system to the desired one
 * The desired endianness is specified by tobig, while
 * the endianness of the current operating system is automatically computed.
 * A conversion takes place only if the two differ.
 * The enidan-corrected variable is saved in output, however,
 * for in-place conversion, pass NULL to output. */
unsigned char *fix_endian(void *output, void *input, size_t len, int tobig)
{
  size_t i, ir;
  static int sysbig = -1;
  unsigned char *fixed = (unsigned char *) output;
  unsigned char *p     = (unsigned char *) input;

  if (sysbig < 0) { /* initial determine the machine's endianess */
    unsigned int feff = 0xFEFF;
    unsigned char *s  = (unsigned char *)&feff;

    sysbig = (*s == 0xFF) ? 0 : 1;
#ifdef ENDIAN_DBG_
    fprintf(stderr, "The current system is %s-endian.\n", 
        sysbig ? "big" : "little");
#endif
  }

  if (tobig == sysbig) { /* no conversion is required */
    if (fixed == NULL) {
      return p;
    } else {
      for (i = 0; i < len; i++) 
        fixed[i] = p[i];
      return fixed;
    }
  } else {
    if (fixed == NULL) { /* in-place conversion */
      for (i = 0; i < (len/2); i++) {
        char ch;
        ir = len - i - 1;
        ch    = p[i];
        p[i]  = p[ir];
        p[ir] = ch;
      }
    } else { /* out-of-place conversion */
      for (i = 0; i < len; i++)
        fixed[len - i - 1] = p[i];
    }
  }
  return p;
}

#endif

