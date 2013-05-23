/* random graph */
#include "clus.c"

double dis_av = 1.0;
double dis_dev = 0.1;
double mu = 0.51; /* 0.5*dis_av is threshold for dis_dev is 0 */

static float **getdismat(int n)
{
  float **mat;
  int i, j;

  mat = dismat_alloc(n);
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      mat[i][j] = dis_av + dis_dev*(2.*rnd0() - 1.);
    }
  }
  return mat;
}

int main(void)
{
  float **mat;
  int n = 300;
  clsys_t *cls;

  mat = getdismat(n);
  cls = cls_init(mat, NULL, n, mu);
  cls_anneal(cls, 3000*n, CLUS_HEATBATH, 1.0, 10.0);
  printf("mu = %g, %d points, %d clusters, ene = %g\n",
      cls->mu0, cls->np, cls->nc, cls->ene);
  cls_write(cls, "clus.txt", NULL, NULL, 0);
  dismat_free(mat);
  cls_free(cls, 1);

  return 0;
}
