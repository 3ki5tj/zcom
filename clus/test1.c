#include "clus.c"

#define MULTI 4
double frac[MULTI] = {1.0, 0.5, 1.0, 1.0};

int weighted = 0;
int n = 200;
int verbose = 2;
int zalgo = 1;
int method = 0; // CLUS_HEATBATH;
double mu = 0.1;
double bet0 = 1.0, bet1 = 10.0;
int nbet = 10;
int nstmin = 1000;

static int dblcmp(const void *pa, const void *pb)
{
  double a = *(double *) pa, b = *(double *) pb;
  return (a > b) ? 1 : (a < b) ? -1 : 0;
}

/* generate a sample with rho(x) = sin(x)^16 */
static double *gensamp1(int n)
{
  double *arr, x, y;
  int i, ix;

  xnew(arr, n);
  for (i = 0; i < n; i++) {
    for (;;) {
      x = rnd0() * MULTI;
      ix = (int) x;
      y = sin(x * M_PI);
      y = y*y;
      y = y*y;
      y = y*y;
      y = y*y;
      if (rnd0() < y * frac[ix]) break;
    }
    arr[i] = x;
  }
  qsort(arr, n, sizeof(*arr), dblcmp);
  return arr;
}

static double *gensamp2(int n)
{
  double *arr, *wt, x, y;
  int i;

  xnew(arr, 2*n);
  wt = arr + n;
  for (i = 0; i < n; i++) {
    x = rnd0() * MULTI;
    y = sin(x * M_PI);
    arr[i] = x;
    wt[i] = y*y*y*y;
  }
  return arr;
}

static float **getdismat(double *arr, int n)
{
  float **mat;
  int i, j;
  
  mat = dismat_alloc(n);
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      mat[i][j] = fabs(arr[i] - arr[j]);
    }
  }
  return mat;
}

int main(void)
{
  double *arr, *wt = NULL;
  float **mat;
  clsys_t *cls;

  if (weighted) {
    arr = gensamp2(n);
    wt = arr + n;
  } else {
    arr = gensamp1(n);
  }
  mat = getdismat(arr, n);
  cls = cls_init(mat, wt, n, mu);
  if (zalgo)
    cls_zalgo(cls, 100*n*n, method, bet0, bet1, nbet, nstmin, verbose);
  else 
    cls_anneal(cls, 10*n*n, method, 0.001, 0.02);
  cls_write(cls, "clus.txt", NULL, NULL, 0);
  printf("Final, energy %g, %d clusters\n", cls->ene, cls->nc);

  {
    FILE *fp;
    int i;
    xfopen(fp, "x.txt", "w", return -1);
    for (i = 0; i < n; i++)
      fprintf(fp, "%d: %g, %g\n", i, arr[i], wt ? wt[i] : 1.0);
    cls->ene = cls_ene(cls, CLUS_CHECK|CLUS_VERBOSE);
    printf("Energy %g\n", cls->ene);
  }

/*
  for (i = 0; i < cls->nc; i++) {
    clus_t *c = cls->c + i;
    for (j = 0; j < c->cnt; j++) {
      double w;
      id = c->idx[j];
      w =  (weighted) ? arr[id+n] : 1.;
      fprintf(stderr, "%g %d %g\n", arr[id], i, w);
    }
  }
*/
  free(arr);
  cls_free(cls, 1);

  return 0;
}
