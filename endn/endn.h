#ifndef ENDN_H__
#define ENDN_H__

#include <stdio.h>
#include <string.h>

int endn_system(void);
size_t endn_fwrite(void *ptr, size_t size, size_t n, FILE *fp, int endn);
size_t endn_fread(void *ptr, size_t size, size_t n, FILE *fp, int flip);
int endn_rmatch(void *src, const void *ref, size_t size, FILE *fp);
int endn_rmatchi(int *src, int iref, FILE *fp);

#endif

