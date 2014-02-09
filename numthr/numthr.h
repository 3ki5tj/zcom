#include "util.h"
#ifndef NUMTHR_H__
#define NUMTHR_H__




INLINE char *getisprime(char *isp, int n)
{
  int i, j, imax;
  for (i = 1; i < n; i++) isp[i] = (char) (i % 2);
  isp[1] = 0;
  imax = (int) sqrt(i);
  for (i = 3; i < imax; i += 2)
    if (isp[i]) /* remove multiples of i */
      for (j = i*i; j < n; j += 2*i)
        isp[j] = 0;
  return isp;
}



INLINE int gcd(int a, int b)
{
  int c;
  while ( a != 0 ) c = a, a = b % a, b = c;
  return b;
}



INLINE int64_t gcd64(int64_t a, int64_t b)
{
  int64_t c;
  while ( a != 0 ) c = a, a = b % a, b = c;
  return b;
}



/* moduler inverse by the Extended Euclidean algorithm
 * http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers */
INLINE int modinv(int a, int n)
{
  int t = 0, r = n, nt = 1, nr = a, t0, r0, q;
  while (nr != 0) {
    q = r / nr;
    t0 = nt; nt = t - q * nt; t = t0;
    r0 = nr; nr = r - q * nr; r = r0;
  }
  if (r > 1) return -1;
  return (t < 0) ? t + n : t;
}



INLINE int64_t modinv64(int64_t a, int64_t n)
{
  int64_t t = 0, r = n, nt = 1, nr = a, t0, r0, q;
  while (nr != 0) {
    q = r / nr;
    t0 = nt; nt = t - q * nt; t = t0;
    r0 = nr; nr = r - q * nr; r = r0;
  }
  if (r > 1) return -1;
  return (t < 0) ? t + n : t;
}



/* http://en.wikipedia.org/wiki/Binary_exponentiation */
INLINE int modpow(int x, int e, int m)
{
  int y;
  for ( y = 1; e > 0; x = x * x % m, e >>= 1 )
    if ( e & 0x1 ) y = y * x % m;
  return y;
}



INLINE int64_t modpow64(int64_t x, int64_t e, int64_t m)
{
  int64_t y;
  for ( y = 1; e > 0; x = x * x % m, e >>= 1 )
    if ( e & 0x1 ) y = y * x % m;
  return y;
}



/* solve r^2 = n (mod p), return one x, the other solution is p - r
 * http://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
 * p must be an odd prime */
INLINE int modpsqrt(int n, int p)
{
  int p1 = (p - 1)/2, q, s, z, c, r, t, m, t2, i, j, b;
  /* compute the Legendre symbol, if it is -1, no quadratic residue */
  if (modpow(n, p1, p) != 1) return -1;
  for ( s = 1, q = p1; q % 2 == 0; s++ ) q /= 2;
  if ( s == 1 ) return modpow(n, (p + 1) / 4, p);
  /* find a non-quadratic residue z */
  for ( z = 2; z < p && modpow(z, p1, p) != p - 1; z++ ) ;
  r = modpow(n, (q + 1) / 2, p);
  c = modpow(z, q, p);
  t = modpow(n, q, p);
  m = s;
  for ( ; t != 1; r = r * b % p, c = b * b % p, t = t * c % p, m = i ) {
    for ( t2 = t, i = 0; i < m && t2 != 1; i++ ) t2 = t2 * t2 % p;
    for ( b = c, j = 0; j < m - i - 1; j++ ) b = b * b % p;
  }
  return r;
}



/* solve r^2 = n (mod p), return one x, the other solution is p - r
 * http://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
 * p must be an odd prime */
int64_t modpsqrt64(int64_t n, int64_t p)
{
  int64_t p1 = (p - 1)/2, q, s, z, c, r, t, m, t2, i, j, b;
  /* compute the Legendre symbol, if it is -1, no quadratic residue */
  if (modpow64(n, p1, p) != 1) return -1;
  for ( s = 1, q = p1; q % 2 == 0; s++ ) q /= 2;
  if ( s == 1 ) return modpow64(n, (p + 1) / 4, p);
  /* find a non-quadratic residue z */
  for ( z = 2; z < p && modpow64(z, p1, p) != p - 1; z++ ) ;
  c = modpow64(z, q, p);
  r = modpow64(n, (q + 1) / 2, p);
  t = modpow64(n, q, p);
  m = s;
  for ( ; t != 1; r = r * b % p, c = b * b % p, t = t * c % p, m = i ) {
    for ( t2 = t, i = 0; i < m && t2 != 1; i++ ) t2 = t2 * t2 % p;
    for ( b = c, j = 0; j < m - i - 1; j++ ) b = b * b % p;
  }
  return r;
}


#endif

