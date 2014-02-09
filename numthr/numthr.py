from math import *


def getdigits(n, p):
  ''' write n as a base-p number
      the least significant digit goes first
  '''
  arr = []
  while n:
    arr += [n % p, ]
    n /= p
  return arr



def gcd(a, b):
  ''' return the greatest common divisor of a and b
      also from fractions import gcd
  '''
  # make sure a <= b
  while a: a, b = b % a, a
  return b



def lcm(a, b):
  return a / gcd(a, b) * b



def xgcd(a, b):
  ''' solve a * s + b * t = gcd(a, b)
      return r, s, t '''
  s, so = 0, 1
  t, to = 1, 0
  r, ro = b, a
  while r:
    q = ro / r
    ro, r = r, ro - q * r
    so, s = s, so - q * s
    to, t = t, to - q * t
  return ro, so, to



def getprimes(n):
  """ return a list of primes <= n
  http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python
  """
  n += 1
  corr = (n % 6 > 1)
  n += [0, -1, 4, 3, 2, 1][n % 6]
  s = [False] + [True]*(n/3 - 1)
  for i in xrange(int(n**.5)/3 + 1):
    if s[i]:
      p = 3*i+1 | 1 # if i is odd, hence 3*i+1 is even, p = 3*i+2, otherwise p = 3*i+1
      s[p*p/3                 : : 2*p] = [False]*((n/6-(p*p)/6-1)/p+1)
      s[(p*p+4*p-2*p*(i&1))/3 : : 2*p] = [False]*((n/6-(p*p+4*p-2*p*(i&1))/6-1)/p+1)
  return [2, 3] + [3*i+1 | 1 for i in xrange(1, n/3 - corr) if s[i]]



def getisprime(n):
  """ return an boolean array to test if numbers are primes
      isp = getisprime()
      if isp[n]: ...
  """
  if n % 2: n += 1
  s = [False, False] + [True]*(n - 1)
  s[4 : : 2] = [False]*(n/2 - 1) # remove multiples of 2
  for p in xrange(3, int(n**.5) + 1, 2):
    # numbers smaller than p*p have been handled by smaller primes
    # further since p is odd, p*p + odd * p is even, and has been removed
    # we only need to count the number of multiples from p*p to n
    if s[p]: s[p*p : : 2*p] = [False]*((n - p*p - 1)/(2*p) + 1)
  return s

#n = 1000
#print getprimes(n),"\n",[x for x in xrange(n) if getisprime(n)[x]]
#
# test timing
#n = 100000000
#import time
#a = time.time(); print len(getprimes(n)); b = time.time(); print "getprimes", b - a
#u = getisprime(n); print "getisprime", time.time() - b



def factors(n, primes = None):
  ''' factors n = p1^e1 p2^e2 ...
      return a dictionary {p1:e1, p2:e2, ..., }
      primes should contain all primes <= sqrt(n)
  '''
  import math
  if not primes: primes = getprimes(int(math.sqrt(n)) + 1)
  dic = {}
  for p in primes:
    if p*p > n: break
    e = 0
    while n % p == 0:
      n /= p
      e += 1
    if e: dic[p] = e
  if n > 1: dic[n] = 1
  return dic

#print factors(10**9 + 1)



def moebiusmu(n, primes = None):
  ''' return Moebius function
        mu(n) = 0 if n contains a prime square, or
              = (-1)^(number of distinct prime factors)
      http://en.wikipedia.org/wiki/Euler_totient_function#Euler.27s_product_formula
      the optional parameter `primes' should contain primes <= n
  '''
  facs = factors(n, primes)
  mu = 1
  for p in facs:
    if facs[p] > 1: return 0
    mu *= -1
  return mu



def eulerphi(n, primes = None):
  ''' return Euler's totient function
      phi(n) = n prod_{p|n} (1 - 1/p) = prod_{p|n} (p^e - p^(e-1))
      http://en.wikipedia.org/wiki/Euler_totient_function#Euler.27s_product_formula
      the optional parameter `primes' should contain primes <= n
  '''
  facs = factors(n, primes)
  phi = 1
  for p in facs:
    e = facs[p]
    phi *= p**(e - 1) * (p - 1)
  return phi



def partialphi(n, primes):
  ''' return the number of integers <= n
      that are not multiples of primes in `primes'
      n - [n/p1] - [n/p2] - ... + [n/p1/p2] + ... '''
  if n <= 0: return 0
  elif n == 1: return 1
  m = len(primes)
  if m == 1: return n - n/primes[0]
  elif m == 2: return n - n/primes[0] - n/primes[1] + n/primes[0]/primes[1]
  elif m == 3:
    q0, p1, p2 = n/primes[0], primes[1], primes[2]
    return n - q0 - n/p1 - n/p2 + q0/p1 + n/p1/p2 + q0/p2 - q0/p1/p2
  elif m == 4:
    q0, p1, p2, p3 = n/primes[0], primes[1], primes[2], primes[3]
    return n - q0 - n/p1 - n/p2 - n/p3 \
        + q0/p1 + n/p1/p2 + q0/p2 + q0/p3 + n/p1/p3 + n/p2/p3 \
        - q0/p1/p2 - q0/p1/p3 - q0/p2/p3 - n/p1/p2/p3 + q0/p1/p2/p3
  sm = 0
  for x in range(1 << m):
    ps = [primes[i] for i in range(m) if x & 1<<i]
    prod = reduce(lambda a, b: a*b, ps, 1)
    sm += (n / prod) * (1 - len(ps) % 2 * 2)
  return sm

#print partialphi(16, [2, 5])



def factorial(n):
  ''' return n! '''
  return reduce(lambda a, b: a*b, range(2, n + 1))


def binomial(n, m):
  if 2*m > n: m = n - m
  x = 1
  for k in range(m):
    x = x * (n - k) / (k + 1)
  return x


def modinv0(a, m):
  ''' return the modular multiplicative inverse, direct method '''
  for x in range(1, m):
    if a * x % m == 1:
      return x
  return 0 # error



def modinv(a, m):
  ''' return the modular multiplicative inverse
      by the extended Euclidean algorithm
      http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Extended_Euclidean_algorithm
      a * x = 1 (mod m) means a * x - m * q = 1
      the last problem can be solved by the extended Euclidean algorithm
      given that gcd(a, m) = 1, otherwise the inverse doesn't exist
      http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
  '''
  #return -xgcd(a % m, -m)[1]
  t, nt, r, nr = 0, 1, m, a
  while nr:
    q = r / nr
    t, nt, r, nr = nt, t - q * nt, nr, r - q * nr
  if r > 1: # not invertible
    raise Exception
  return t if t > 0 else t + m

#print modinv0(3, 4), modinv(3, 4)


def modinv_euler(a, m, primes = None):
  ''' return the modular multiplicative inverse using Euler's theorem
      http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Using_Euler.27s_theorem
      since a ^ phi(m) = 1 (mod m), a^(-1) = a^(phi(m) - 1) (mod m)
      here phi(m) is Euler's totient function
      the optional parameter primes should contain primes <= a
  '''
  # return modpow(a, eulerphi(a) - 1, m)
  # the system function pow() takes the last argument to mean modulo m
  return pow(a, eulerphi(a) - 1, m)



def modpow0(x, e, m):
  ''' return x^e mod m, direct method '''
  y = 1
  for i in range(e):
    y = y * x % m
  return y



def modpow(x, e, m):
  ''' return x^e mod m, right-to-left binary method
      http://en.wikipedia.org/wiki/Modular_exponentiation#Right-to-left_binary_method
  '''
  y = 1
  while e > 0:
    if e % 2:
      y = y * x % m
    e = e >> 1 # e = e // 2
    x = x * x % m
  return y


def chineseremainder2(a, m, b, n):
  ''' if x = a (mod m) = b (mod n) with gcd(m, n) = 1
      return x mod m*n from the solution
      x = a * n * (invn)_m + b * m * (invm)_n
      http://en.wikipedia.org/wiki/Chinese_remainder_theorem#Case_of_two_equations
  '''
  invm = modinv(m, n)
  invn = modinv(n, m)
  print m, n, invm, invn
  return (a * n * invn + b * m * invm) % (m * n)

#print chineseremainder2(2, 3, 3, 4)


def chineseremainder(a, m):
  ''' array case of the above
      http://en.wikipedia.org/wiki/Chinese_remainder_theorem#General_case
  '''
  M = 1
  for mi in m: M *= mi
  x = 0
  for i in range(len(a)):
    ai = a[i]
    Mi = M / m[i]
    x = (x + ai * Mi * modinv(Mi, m[i])) % M
  return x

# the classic example in Mathematical Treatise in Nine Sections
#print chineseremainder([2, 3, 2], [3, 5, 7]), "should be", 23
# the wikipedia example
#print chineseremainder([2, 3, 1], [3, 4, 5]), "should be", 11

