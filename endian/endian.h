#ifndef ENDIAN_H__
#define ENDIAN_H__

#include <stdio.h>

unsigned char *fix_endian(void *output, void *input, size_t len, int tobig);

#endif

