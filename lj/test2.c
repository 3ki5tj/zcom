#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 256;
double rho = 0.8, rc = 2.5, L, *x, *v, *f, fs[900];

/* initialize a fcc lattice */
static void initfcc(void)
{
  int i, j, k, id, n1;
  double a;

  n1 = (int) (pow(2*N, 1.0/3) + .999999); /* # of particles per side */
  a = 1./n1;
  for (id = 0, i = 0; i < n1 && id < N; i++)
    for (j = 0; j < n1 && id < N; j++)
      for (k = 0; k < n1 && id < N; k++) {
        if ((i+j+k) % 2 != 0) continue;
        x[id*3 + 0] = (i+.5)*a;
        x[id*3 + 1] = (j+.5)*a;
        x[id*3 + 2] = (k+.5)*a;
        id++;
      }
}

/* compute force */
static void myforce(double *f)
{
  double dx[3], dr2, dr6, fs, rc2 = rc*rc;
  int i, j, d;

  for (i = 0; i < N*3; i++) f[i] = 0;
  for (i = 0; i < N - 1; i++) {
    for (j = i+1; j < N; j++) {
      dx[0]=(x[i*3] - x[j*3])*L;
      dx[1]=(x[i*3+1] - x[j*3+1])*L;
      dx[2]=(x[i*3+2] - x[j*3+2])*L;
      dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; /* no pbc for simplicity */
      if (dr2 > rc2) continue;
      dr2 = 1.f/dr2;
      dr6 = dr2*dr2*dr2;
      fs = dr6*(48.f*dr6-24.f)*dr2;
      for (d = 0; d < 3; d++) {
        f[i*3+d] += dx[d]*fs;
        f[j*3+d] -= dx[d]*fs;
      }
    }
  }
}

#define xnew(a, n) a = calloc(n, sizeof(*(a)))
int main(void)
{
  int t;

  xnew(x, N*3);
  xnew(v, N*3);
  xnew(f, N*3);
  L = pow(1.*N/rho, 1./3);
  initfcc();
  for (t = 1; t <= 40000; t++) {
#ifdef STACK
    myforce(fs);
#else
    myforce(f);
#endif
  }
  return 0;
}
