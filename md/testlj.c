#include "md.h"
#include "include/av.h"

#define D        3
#define N        256
#define DOF      (D*N-D*(D+1)/2)  /* degrees of freedom */
double rho = 0.8;
double T = 1.0;
double rc = 2.5;
double dt = 0.002; /* time step for Newton's equation */
double L;
double x[N][3]; /* from 0 to 1 */
double v[N][3], f[N][3]; /* velocity and force */
double U, Us, K, Vir, p;
double ushift; /* potential energy shift due to truncation, used to check energy conservation */
double Utail; /* tail correction for potential energy */
double ptail; /* tail correction for pressure */
int thermostat = 1; /* use a thermostat */
double thermdt = 0.01; /* time step for thermostat */
int step = 0;
av_t avU, avK, avp;

/*  to the nearest image */
INLINE double pbc(double a)
  { return L * (a - ((int)((a)+1000.5) - 1000)); }

INLINE double *vpbc(double a[])
  { a[0] = pbc(a[0]); a[1] = pbc(a[1]); a[2] = pbc(a[2]); return a; }

INLINE double pbcdist2(double dx[], const double a[], const double b[])
  { return rv3_sqr(vpbc(rv3_diff(dx, a, b))); }

/* initialize an md system */
static void initmd(void)
{
  int i, j, k, n1, id;
  double a, irc3, irc6;

  L = pow(N/rho, 1.0/D);
  if ((rc = 2.5) > L*.5) rc = L*.5;
  irc3 = 1.0/(rc*rc*rc); irc6 = irc3*irc3;
  ushift = 4*irc6*(irc6-1);
  Utail = 8*M_PI*rho*N/9*(irc6 - 3)*irc3;
  ptail = 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3;
  printf("T %.2f, rho %.2f, N %d, L %6.3f, ushift %6.3f, utail %6.3f, ptail %6.3f\n",
    T, rho, N, L, ushift, Utail/N, ptail);

  /* init. a fcc lattic */
  n1 = (int) (pow(2*N, 1.0/D) + .999999); /* # of particles per side */
  a = 1./n1;
  for (id = 0, i = 0; i < n1 && id < N; i++)
    for (j = 0; j < n1 && id < N; j++)
      for (k = 0; k < n1 && id < N; k++) {
        if ((i+j+k) % 2 != 0) continue;
        x[id][0] = (i+.5)*a;
        x[id][1] = (j+.5)*a;
        x[id][2] = (k+.5)*a;
        id++;
      }

  /* init. random velocities */
  for (id = 0; id < N; id++)
    for (j = 0; j < D; j++) v[id][j] = rnd0() - .5;

  md_shiftcom3d(x, N);
  md_shiftang3d(x, v, N);
}

/* compute force and virial, return energy */
static void force(void)
{
  double dx[3], fi[3], dr2, dr6, fs, tmp;
  int i, j, d, prcnt = 0;

  for (i = 0; i < N; i++) f[i][0] = f[i][1] = f[i][2] = 0;
  for (U = Vir = 0, i = 0; i < N - 1; i++) {
    for (d = 0; d < D; d++) fi[d] = 0;
    for (j = i+1; j < N; j++) {
      dr2 = pbcdist2(dx, x[i], x[j]);
      if (dr2 > rc*rc) continue;
      dr2 = 1/dr2;
      dr6 = dr2*dr2*dr2;
      fs = dr6*(48*dr6-24); /* f.r */
      Vir += fs;
      fs *= dr2;
      for (d = 0; d < D; d++) {
        tmp = dx[d]*fs;
        fi[d] += tmp;
        f[j][d] -= tmp;
      }
      U += 4*dr6*(dr6-1);
      prcnt++;
    }
    for (d = 0; d < D; d++)
      f[i][d] += fi[d];
  }
  Us = U - prcnt*ushift; /* shifted energy */
  U += Utail; /* unshifted energy */
  p = rho*T+Vir/(3*L*L*L)+ptail;
}

/* velocity verlet */
static void vv(void)
{
  int i;
  double dth = dt*.5, dtl = dt/L;

  for (i = 0; i < N; i++) { /* VV part 1 */
    rv3_sinc(v[i], f[i], dth);
    rv3_sinc(x[i], v[i], dtl);
  }
  force(); /* calculate the new force */
  for (i = 0; i < N; i++) /* VV part 2 */
    rv3_sinc(v[i], f[i], dth);

  K  = md_ekin3d(v, N, DOF, NULL);
  if (thermostat) md_vrescale3d(v, N, DOF, T, thermdt, &K, NULL);
  step++;
}

int main(void)
{
  int i, nmax = 20000;
  double u, k;

  initmd();
  for (i = 0; i < nmax; i++) {
    vv();
    if (i > nmax/2) { av_add(&avU, U); av_add(&avK, K); av_add(&avp, p); }
  }
  u = av_getave(&avU)/N;
  k = av_getave(&avK)/N;
  printf("U/N = %6.3f, K/N = %6.3f, p = %6.3f\n", u, k, av_getave(&avp));
  return 0;
}

