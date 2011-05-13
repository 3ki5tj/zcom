#include <stdlib.h>
#include "savgol.c"

static void test1d(int w, int ord)
{
  int i;
  double *c;
  FILE *fp = fopen("savgol1d.dat", "w");
  c = savgol(w, ord, 0, 0, 0);
  for (i = 0; i <= 2*w; i++)
    fprintf(fp, "%d %g\n", i-w, c[i]);
  free(c);
  fclose(fp);
}

static void test2d(int w, int ord)
{
  int i, j, id;
  double *c, sum = 0.;
  FILE *fp = fopen("savgol2d.dat", "w");
  c = savgol2d(w, w, ord, 0, 0);
  for (i = 0; i <= 2*w; i++) {
    for (j = 0; j <= 2*w; j++) {
      id = i*(2*w+1)+j;
      sum += c[id];
      fprintf(fp, "%d %d %g\n", i - w, j - w , c[id]);
    }
    fprintf(fp, "\n");
  }
  printf("2d savgol sum = %g\n", sum);
  free(c);
  fclose(fp);
}

int main(int argc, char **argv)
{
  double *c, *c2;
  int w = 2, ord = 2;

  if (argc > 1) w = atoi(argv[1]);
  c = savgol(w, ord, 0, 0, 1);
  printf("\n\n");
  c2 = savgol2d(w, w, ord, 0, 1);
  printf("\n");
  free(c);
  free(c2);

  test1d(4000, ord);
  test2d(20, ord);
  return 0;
}
