#ifndef ENDN_H__
#define ENDN_H__

#include <stdio.h>
#include <string.h>

int endn_system(void);

#define endn_converti(p, size, n, endn) endn_convert(NULL, p, size, n, endn)
void *endn_convert(void *dest, void *src, size_t size, size_t n, int endtar);

#endif

