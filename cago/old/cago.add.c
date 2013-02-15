/* the change of the potential energy */
INLINE real cago_depot0(cago_t *go, rv3_t *x, int i, rv3_t xi)
{
  /* stupid implementation, assuming go->epot is correct */
  rv3_t xo;
  real ep0 = go->epot, ep1;

  rv3_copy(xo, x[i]);
  rv3_copy(x[i], xi);
  ep1 = cago_force(go, x, NULL);
  rv3_copy(x[i], xo); /* recover the displaced particle */
  return ep1 - ep0;
}

/* metropolis algorithm (with periodic checking) */
INLINE int cago_metro1(cago_t *go, real amp, real bet)
{
  int i, j;
  rv3_t xi;
  real du;
  static int cc;

  i = (int) (go->n * rnd0());
  for (j = 0; j < 3; j++) {
    xi[j] = go->x[i][j] + (rnd0() * 2.f - 1) * amp;
  }
  du = cago_depot(go, go->x, i, xi);
  if (du < 0 || rnd0() < exp(-bet * du)) {
    rv3_copy(go->x[i], xi);
    go->epot += du;
    if (++cc % 1000 == 0) { /* refresh energy */
      real ep = cago_force(go, go->x, NULL);
      if (fabs(ep - go->epot) > 0.01) {
        fprintf(stderr, "energy corruptions %g vs. %g\n", go->epot, ep);
      }
      go->epot = ep;
    }
    return 1;
  } else {
    return 0;
  }
}


