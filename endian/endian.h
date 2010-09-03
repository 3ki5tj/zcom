#ifndef ENDIAN_H__
#define ENDIAN_H__

#include <stdio.h>
#include <string.h>

#define fix_endian_inp(p, size, n, endn) fix_endian(NULL, p, size, n, endn)
void *fix_endian(void *dest, void *src, size_t size, size_t n, int endtar);

#endif

