#include "rng.h"
#include "util.c"

#ifndef POTTS2_C__
#define POTTS2_C__

#include "potts2.h"

/* compute the total energy and magnetization */
int pt2_em(potts_t *pt)
{
  int i, j, l, s, s1, s2, *p;

  pt->E = 0;
  p = pt->s;
  l = pt->l;
  for (i = 0; i < pt->q; i++) pt->M[i] = 0;
  for (i = 0; i < l; i++)
    for (j = 0; j < l; j++) {
      s = p[i*l + j];
      s1 = p[((i+1)%l)*l + j];
      if (s1 == s) pt->E--;
      s2 = p[i*l + (j+1)%l];
      if (s2 == s) pt->E--;
      pt->M[s]++;
    }
  return pt->E;
}

int pt2_check(potts_t *pt)
{
  int i, e, *mg, q;

  q = pt->q;
  for (i = 0; i < pt->n; i++) /* check spin value */
    if (pt->s[i] < 0 || pt->s[i] >= q) {
      fprintf(stderr, "error: s[%d] = %d\n", i, pt->s[i]);
      return -1;
    }
  e = pt->E;
  xnew(mg, pt->q);
  for (i = 0; i < pt->q; i++)
    mg[i] = pt->M[i];
  if (e != pt2_em(pt)) { /* check energy */
    fprintf(stderr, "error: E = %d, should be %d\n", 
        e, pt->E);
    free(mg);
    return -1;
  }
  for (i = 0; i < pt->q; i++) {
    if (mg[i] != pt->M[i]) {
      fprintf(stderr, "error: M(%d) = %d, should be %d",
          i, mg[i], pt->M[i]);
      free(mg);
      return -1;
    }
  }
  free(mg);
  return 0;
}

/* pick a random site (return its id)
 * compute h[j], the numbers of neighboring spins with value j */
/* load spin configuration */
int pt2_load(potts_t *pt, const char *fname)
{
  FILE *fp;
  int i, lx, ly, n, c;
  char s[80];

  if ((fp = fopen(fname, "r")) == NULL) {
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "missing first line %s\n", fname);
    return -1;
  }
  if (4 != sscanf(s, "%d%d%d%d", &i, &lx, &ly, &n) 
      || i != 2 || lx != ly || lx != pt->l || n != pt->n) {
    fprintf(stderr, "bad setting: %dD, %dx%d = %d\n", i, lx, ly, n);
    return -1;
  }
  for (i = 0; i < n; i++) {
    while ((c=fgetc(fp)) != EOF && c == '\n') ;
    if (c == EOF) break;
    c -= '0';
    if (c < 0 || c >= pt->q) {
      fprintf(stderr, "BAD %s s[%d] = %d, q = %d\n", fname, i, c, pt->q);
      break;
    }
    pt->s[i] = c;
  }
  if (i < n) {
    fprintf(stderr, "%s: data stopped at i = %d, clear\n", fname, i);
    for (i = 0; i < n; i++) pt->s[i] = 0;
  }
  fclose(fp);
  pt2_em(pt); /* re-compute energy/magnetization */
  return 0;
}

/* save spin configuration */
int pt2_save(const potts_t *pt, const char *fname)
{
  FILE *fp;
  int i, j, l, *p;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  l = pt->l;
  fprintf(fp, "%d %d %d %d\n", pt->d, l, l, pt->n);
  for (p = pt->s, i = 0; i < l; i++) {
    for (j = 0; j < l; j++, p++)
      fprintf(fp, "%c", *p+'0');
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* initialize an lxl q-state Potts model */
potts_t *pt2_open(int l, int q)
{
  int i, n;
  potts_t *pt;

  xnew(pt, 1);
  pt->d = 2;
  pt->q = q;
  pt->l = l;
  pt->n = n = l*l;
  xnew(pt->s, n);
  xnew(pt->M, q);
  for (i = 0; i < n; i++)
    pt->s[i] = 0;
  for (pt->M[0] = n, i = 1; i < q; i++)
    pt->M[i] = 0;
  pt->E = -pt->d * n;
  xnew(pt->accprb, q+1);
  pt->accprb[0] = 0.;
  /* dynamic array of uproba/dproba seems to be faster */
  xnew(pt->uproba, 2*pt->d+1);
  pt->uproba[0] = 0xffffffffu;
  xnew(pt->dproba, 2*pt->d+1);
  pt->dproba[0] = 1.;
  return pt;
}

void pt2_close(potts_t *pt)
{
  if (pt != NULL) {
    free(pt->s);
    free(pt->M);
    free(pt->accprb);
    free(pt->uproba);
    free(pt->dproba);
    free(pt);
  }
}

#endif

