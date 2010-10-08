#include "rng.h"
#ifndef IS2_C__
#define IS2_C__

#include "is2.h"

/* compute total energy and magnetization */
int is2_em(is_t *is)
{
  int l, i, j, s, u, e, m, *p, *pu;

  e = m = 0;
  p = is->s;
  l = is->l;
  for (i = 0; i < l; ) {
    pu = (++i == l) ? is->s : p+l;
    for (j = 0; j < l; ) {
      m += (s = *p++);
      u = *pu++;
      e += s*(u + ((++j == l) ? *(p-l) : *p));
    }
  }
  is->M = m;
  return is->E = -e;
}

int is2_check(is_t *is)
{
  int i, e, m;

  for (i = 0; i < is->n; i++) 
    if (is->s[i] != 1 && is->s[i] != -1) {
      fprintf(stderr, "error: s[%d] = %d\n", i, is->s[i]);
      return -1;
    }
  e = is->E;
  m = is->M;
  is2_em(is);
  if  (e != is->E || e < -2*is->n || e > 2*is->n
    || m != is->M || m < -is->n   || m > is->n) {
    fprintf(stderr, "error: E = %d, %d; M = %d, %d\n", 
        e, is->E, m, is->M);
    return -1;
  }
  return 0;
}

/* pick a random site, count neighbors with different spins */
int is2_pick(const is_t *is, int *h)
{
  int id, ix, iy, l, lm, n, nm, *p;

  lm = (l = is->l) - 1;
  nm = (n = is->n) - l;
  id = (int)(rnd0() * n);
  iy = id / l, ix = id % l;
  p = is->s + id;
  *h = *p * ( ((ix != 0 ) ? *(p-1) : *(p+lm))   /* left  */
            + ((ix != lm) ? *(p+1) : *(p-lm))   /* right */
            + ((iy != 0 ) ? *(p-l) : *(p+nm))   /* down  */
            + ((iy != lm) ? *(p+l) : *(p-nm))); /* up    */
  return id;
}

/* flip site id, with nd different neighbors */
int is2_flip(is_t *is, int id, int h)
{
  assert(id < is->n);
  is->M += (is->s[id] = -is->s[id])*2;
  return is->E += h*2;
}

int is2_load(is_t *is, const char *fname)
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
      || i != 2 || lx != ly || lx != is->l || n != is->n) {
    fprintf(stderr, "bad setting: %dD, %dx%d = %d\n", i, lx, ly, n);
    return -1;
  }
  for (i = 0; i < n; i++) {
    while ((c=fgetc(fp)) != EOF && c == '\n') ;
    if (c == EOF) break;
    is->s[i] = (c == ' ') ? -1 : 1;
  }
  if (i < n)
    fprintf(stderr, "%s: data stopped at i = %d\n", fname, i);
  fclose(fp);
  is2_em(is);
  return 0;
}

int is2_save(const is_t *is, const char *fname)
{
  FILE *fp;
  int i, j, l, *p;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  l = is->l;
  fprintf(fp, "%d %d %d %d\n", is->d, l, l, is->n);
  for (p = is->s, i = 0; i < l; i++) {
    for (j = 0; j < l; j++, p++)
      fprintf(fp, "%c", (*p > 0) ? '#' : ' ');
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* initialize an lxl Ising model */
is_t *is2_open(int l)
{
  int i, n;
  is_t *is;

  if ((is = calloc(1, sizeof(*is))) == NULL){
    fprintf(stderr, "no memory for is.\n");
    return NULL;
  }
  is->d = 2;
  is->l = l;
  is->n = n = l*l;
  if ((is->s = malloc(sizeof(is->s[0])*n)) == NULL) {
    fprintf(stderr, "no memory for spin, %dx%d\n", l, l);
    return NULL;
  }
  for (i = 0; i < n; i++) is->s[i] = -1;
  is->M = -n;
  is->E = -2*n;
  return is;
}

void is2_close(is_t *is) 
{
  if (is != NULL) {
    free(is->s);
    free(is);
  }
}

#define LN_BIG 50.0
/* c = log(exp(a) + exp(b)) */
static double lnadd_(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a-b) > LN_BIG) ? a : a + log(1 + exp(-c));
}
/* c = log(exp(a)-exp(b)), only works for a>b */
static double lnmin_(double a, double b)
{
  double c;
  assert(a >= b);
  return ((c = a-b) > LN_BIG) ? a : a + log(1 - exp(-c));
}
/* c = log(exp(a)+b) */
static double lnaddn_(double a, double b)
{
  return (a > LN_BIG) ? a : a + log(1 + b*exp(-a));
}
/* exact solution of ising model */
double is2_exact(is_t *is, double beta, double *eav, double *cv)
{
  double lxh, n, ex, f, th, sech, bet2, bsqr, log2, x;
  double lnz, lnz1, lnz2, lnz3, lnz4, dz, ddz;
  double z21, z31, z41, za1;
  double dr1, dr2, dr3, dr4, ddr1, ddr2, ddr3, ddr4;
  double g, g0, dg, ddg, dg0;
  double xn2b, sh2b, coth2b;
  double lnch2b, lncc2b, lncl, lnsl, cd, cdsqr, lnddcl;
  int r, sgn4 = 1, lx, ly;

  lx = is->l;
  ly = is->l;
  lxh = .5*lx;
  n = lx * ly;
  log2 = log(2.0);
  bet2 = 2.*beta;
  bsqr = beta*beta;
  xn2b = exp(-bet2);
  if (lx == 2 && ly == 2) { /* 2x2 system */
    double lnc, lnd;
    x = 8.*beta;
    lnc = lnadd_(x, -x);
    lnd = lnaddn_(lnc, 6.);
    lnz = lnd + log2;
    *eav = -8.*exp(lnmin_(x, -x) - lnd); /* -8*sinh(8*b)/(3+cosh(8*h)) */
    *cv = bsqr * 384. * exp(lnaddn_(lnc,2./3) - 2.0*lnd); /* 64*(1+3cosh(8*b))/(3+cosh(8*b))^2 */
    return lnz;
  } else if (fabs(beta) < 1e-6) { /* high T approx. normal branch unstable if beta < 1e-6 */
    lnz = n * (2.*lnadd_(beta, -beta) - log2);
    x = 1. + xn2b;
    *eav = -2. * n * (1. - xn2b)/x;
    *cv = bsqr * 8.*n*xn2b/(x*x);
    return lnz; /* +n*tanh(beta)^4 */
  }

  lnz1 = lnz2 = lnz3 = lnz4 = 0;
  dr1 = dr2 = dr3 = dr4 = 0;
  ddr1 = ddr2 = ddr3 = ddr4 = 0;
  lnch2b = lnadd_(bet2, -bet2) - log2;
  coth2b = 2./(1. - xn2b*xn2b) - 1.;
  lncc2b = lnch2b + log(coth2b); /* ln[ cosh(2b) * coth(2b) ] */
  g0 = bet2 + log(2./(1. + xn2b) - 1.);
  sgn4 = (g0 >= 0) ? 1 : -1;

  sh2b = 0.5*(1./xn2b - xn2b);
  dg0 = 2. + 2./sh2b;
  x = sh2b*sh2b;
  cd = 2. - 2./x; /* cl' = cd * cosh(2b) */
  cdsqr = cd*cd;
  lnddcl = lnaddn_(lncc2b, 2.0/(x * sh2b)) + 2.*log2; /* log(cl'') */

  for (r = 0; r < ly; r++) { /* for odd number */
    lncl = lnaddn_(lncc2b, -cos((2.*r + 1.)*M_PI/ly));
    lnsl = lncl + 0.5*log(1. - exp(-2.*lncl));
    g = lnadd_(lncl, lnsl);
    f = lxh*g;
    lnz1 += lnadd_(f, -f);
    lnz2 += lnmin_(f, -f);

    dg = exp(lnch2b - lnsl)*cd; /* g' = cl'/sl; */
    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    x = lxh*dg;
    dr1 += x*th;
    dr2 += x/th;

    /* g''=cl''/sl - cl' ^2 *cl/sl^3; */
    ddg = exp(lnddcl - lnsl);
    ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    sech = 2.0*dg/(ex + 1.0/ex); /* g' * sech(0.5*lx*g) */
    ddr1 += lxh*(ddg*th + lxh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex); /* g' * sech(0.5*lx*g) */
    ddr2 += lxh*(ddg/th - lxh*(sech*sech));

    if (r == 0) {
      g = g0;
    } else {
      lncl = lnaddn_(lncc2b, -cos(2.0*M_PI*r/ly));
      lnsl = lncl+0.5*log(1-exp(-2*lncl));
      g = lnadd_(lncl, lnsl);
      assert(g > 0);
    }
    f = lxh*g;
    lnz3 += lnadd_(f, -f); /* log [2 cosh(f)] */
    lnz4 += (f < 0) ? lnmin_(-f, f) : lnmin_(f, -f); /* avoid neg. g0 */
   
    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    dg = (r == 0) ? dg0 : exp(lnch2b - lnsl)*cd;
    dr3 += lxh*dg*th;
    dr4 += lxh*dg/th;

    if (r == 0) {
      ddg = -4*coth2b*coth2b*exp(-lnch2b);
    } else {
      ddg = exp(lnddcl - lnsl);
      ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    }
    sech = 2.0*dg/(ex + 1.0/ex);
    ddr3 += lxh*(ddg*th + lxh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex);
    ddr4 += lxh*(ddg/th - lxh*(sech*sech));
  }

  z21 = exp(lnz2 - lnz1);
  z31 = exp(lnz3 - lnz1);
  z41 = sgn4*exp(lnz4 - lnz1);
  za1 = 1.0 + z21 + z31 + z41;
  lnz = lnz1 + log(za1);
  lnz += .5*n*log(2.*sh2b) - log2;
  dz = (dr1 + z21*dr2 + z31*dr3 + z41*dr4)/za1;
  *eav = - n*coth2b - dz;
  ddr1 += dr1*dr1;
  ddr2 += dr2*dr2;
  ddr3 += dr3*dr3;
  ddr4 += dr4*dr4;
  ddz = (ddr1 + z21*ddr2 + z31*ddr3 + z41*ddr4)/za1;
  *cv = bsqr * (-2.*n/(sh2b*sh2b) + ddz - dz*dz);
  return lnz;
}

#endif /* IS2_C__ */

