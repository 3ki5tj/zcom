/* add spin j to the queue if s[j] is different from s
 * return the spin */
INLINE void is2_addtoqueue_ne(ising_t *is, int j, int s,
    double r, int *cnt)
{
  if ( is->s[j] == s && rand01() < r ) {
    is->queue[ (*cnt)++ ] = j;
    is->s[j] = -s;
  }
}



/* Wolff algorithm (without energy update) */
INLINE int is2_wolff_ne(ising_t *is, double padd)
{
  int l = is->l, n = is->n, i, ix, iy, id, s, cnt = 0;

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  s = is->s[id];
  is->s[id] = -s;
  is->queue[ cnt++ ] = id;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = is->queue[i];
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    is2_addtoqueue_ne(is, iy + (ix + 1) % l,     s, padd, &cnt);
    is2_addtoqueue_ne(is, iy + (ix + l - 1) % l, s, padd, &cnt);
    is2_addtoqueue_ne(is, (iy + l) % n + ix,     s, padd, &cnt);
    is2_addtoqueue_ne(is, (iy + n - l) % n + ix, s, padd, &cnt);
  }

  is->M -= 2 * s * cnt;
  return 0;
}



static void runwolff_ne(ising_t *is, double steps, double beta, int ncheck)
{
  double t, acc, tot, padd;
  double s1, e, se, se2, eav, cv;
  double eref, cvref, lnzref;
  int nt = ncheck;

  padd = 1 - exp(-2*beta);
  acc = tot = 1e-8;
  se = se2 = 0.0;
  for (t = 1.0; t <= steps; t += 1.0) {
    is2_wolff_ne(is, padd);
    if ( --nt == 0 ) {
      is2_em(is);
      s1 += 1;
      se += e = is->E;
      se2 += e * e;
      nt = ncheck;
    }
  }
  eav = se / s1;
  cv = (beta*beta) * (se2/s1 - eav*eav);
  lnzref = is2_exact(is, beta, &eref, &cvref);
  printf("ar: %g, eav: %.6f (%.6f), cv: %.3f (%.3f), lnz: %.6f\n",
      acc/tot, eav, eref, cv, cvref, lnzref);
}





