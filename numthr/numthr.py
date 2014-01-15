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



def factors(n, primes = None):
  ''' factors n = p1^e1 p2^e2 ...
      return a dictionary {p1:e1, p2:e2, ..., }
      primes should contain all primes <= n
  '''
  if not primes: primes = getprimes(q)
  dic = {}
  for p in primes:
    if p > n: break
    e = 0
    while n % p == 0:
      n /= p
      e += 1
    if e: dic[p] = e
  return dic



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
  '''
  return xgcd(a, -m)[1]
  



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
