#include "hist.c"
#include "test/rv3.h"

#define XMIN 0.0
#define XMAX 2.0
#define XDEL 0.01

#define N 1000
double x[N][3];

typedef struct {
  int nfr;
  int n;
  double vol;
} den_t;

static void genhist(hist_t *hs, int nfr, double l, den_t *den)
{
  int i, j, k, ifr;
  double dx[3], r;

  for (ifr = 0; ifr < nfr; ifr++) {
    for (i = 0; i < N; i++) {
      x[i][0] = l*rand()/RAND_MAX;
      x[i][1] = l*rand()/RAND_MAX;
      x[i][2] = l*rand()/RAND_MAX;
    }
    
    for (i = 0; i < N; i++)
      for (j = i+1; j < N; j++) {
        rv3_diff(dx, x[i], x[j]);
        for (k = 0; k < 3; k++) {
          if (dx[k] > .5*l) dx[k] -= l;
          else if (dx[k] < -.5*l) dx[k] += l;
        }
        r = rv3_norm(dx);
        hs_add(hs, &r, 1., HIST_VERBOSE);
      }
  }
  den->nfr = nfr;
  den->n = N;
  den->vol = l*l*l;
}

static int fwheader(FILE *fp, void *pdata)
{
  den_t *den = (den_t *)pdata;  
  fprintf(fp, "RDF %d %d %g | ", den->nfr, den->n, den->vol);
  return 0;
}

static int frheader(const char *s, void *pdata)
{
  den_t *den = (den_t *)pdata;  
  int ret = sscanf(s, " RDF %d%d%lf | ", &(den->nfr), &(den->n), &(den->vol));
  return (ret == 3) ? 0 : 1;
}

static double rdfnorm(int j, int i, double xmin, double dx, void *pdata)
{
  den_t *den = (den_t *)pdata;
  int npr;
  double x, fac, vsph;
 
  (void)j; 
  x = xmin + i*dx;
  vsph = (4.*M_PI/3)*dx*(3*x*(x+dx) + dx*dx);
  npr = den->n*(den->n-1)/2;
  fac = den->vol/vsph/npr/den->nfr;
  return fac;
}

int main(void)
{
  den_t den[1];
  hist_t *hs;

  hs = hs_initx(1, XMIN, XMAX, XDEL, fwheader, frheader, rdfnorm);

  /* generate histogram */
  genhist(hs, 17, 1., den);
  
  hs_savex(hs, "RDF", den, HIST_ADDAHALF|HIST_KEEPHIST);
  hs_savex(hs, "rdf.dat", den, HIST_ADDAHALF);

  /* now try to load histogram */
  if (0 != hs_loadx(hs, "RDF", den, HIST_VERBOSE)) {
    fprintf(stderr, "cannot load histogram\n");
    return -1;
  }
  /* write again */
  hs_savex(hs, "RDF2", den, HIST_ADDAHALF|HIST_KEEPHIST);
  hs_savex(hs, "rdf2.dat", den, HIST_ADDAHALF);
  hs_free(hs);
  return 0; 
}

