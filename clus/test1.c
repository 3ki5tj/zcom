#include "clus.c"

#define MULTI 4
double frac[MULTI] = {1.0, .1, 0.6, 0.3};

static double *gensamp1(int n)
{
  double *arr, x, y;
  int i, ix;

  xnew(arr, n);
  for (i = 0; i < n; i++) {
    for (;;) {
      x = (rnd0()*MULTI)*M_PI;
      ix = (int)(x/M_PI);
      if (rnd0() > frac[ix]) continue;
      y = sin(x);
      y = y*y*y*y;
      if (rnd0() < y) break;
    }
    arr[i] = x;
  }
  return arr;
}

static double *gensamp2(int n)
{
  double *arr, *wt, x, y;
  int i;

  xnew(arr, 2*n);
  wt = arr + n;
  for (i = 0; i < n; i++) {
    x = (rnd0()*MULTI)*M_PI;
    y = sin(x);
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
  int weighted = 0;
  int n = 200;
  clsys_t *cls;
  int verbose = 2;
  int zalgo = 1;
  int method = 0; // CLUS_HEATBATH;
  double mu = 0.2;
  double bet0 = 1., bet1 = 1.;

  if (weighted) {
    arr = gensamp2(n);
    wt = arr + n;
  } else {
    arr = gensamp1(n);
  }
  mat = getdismat(arr, n);
  cls = cls_init(mat, wt, n, mu);
  if (zalgo)
    cls_zalgo(cls, 10*n*n, method, bet0, bet1, 10, verbose);
  else 
    cls_anneal(cls, 1000*n*n, method, 2.0, 10.0);
  cls_write(cls, "clus.txt", NULL, NULL, 0);
  //for (i = 0; i < n; i++)
  //  printf("%d: %g, %g\n", i, arr[i], arr[i+n]);
  printf("energy %g, %d clusters\n", cls->ene, cls->nc);
  /*
  for (i = 0; i < cls->nc; i++) {
    clus_t *c = cls->c + i;
    for (j = 0; j < c->cnt; j++) {
      double w;
      id = c->idx[j];
      w =  (weighted) ? arr[id+n] : 1.;
      fprintf(stderr, "%g %d %g\n", arr[id], i, w);
    }
  }*/
  free(arr);
  cls_free(cls, 1);

  return 0;
}
