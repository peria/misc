#!/usr/bin/python

import fractions

def getPrimes(n):
  isprime = [True for _ in xrange(n / 2)]
  primes = [2]
  for k in xrange(1, n / 2):
    if not isprime[k]:
      continue
    p = k * 2 + 1
    primes.append(p)
    for p2 in xrange(p * p / 2, n / 2, p):
      isprime[p2] = False
  return primes


def factor(primes, n):
  # Use p-1 method
  k, a = 0, 2
  for p in primes:
    if p < 10000:
      a = pow(a, p * 10, n)
      continue
    a = pow(a, p, n)
    k = k + 1
    if k % 100000 != 0:
      continue
    g = fractions.gcd(a - 1, n)
    if 1 < g:
      return g
  return fractions.gcd(a - 1, n)


def test(primes, n):
  p = factor(primes, n)
  if p == n:
    print 'Too many primes to factor %x' % n
  elif 1 == p:
    print 'Fail to factor %x' % n
  else:
    print '%x = %x * %x' % (n, p, n / p)
    print '%d = %d * %d' % (n, p, n / p)


if __name__ == '__main__':
  n = 100000000
  primes = getPrimes(n)
  print 'generated primes up to %d' % n

  # These numbers are in http://www.gsic.titech.ac.jp/supercon/supercon2004/kaisetsu/
  test(primes, 0xd44cdd3dc3436e0c04cdedb0f)  # Contest
  test(primes, 0x851bed3b13e22bd06bc3a0ba7)  # Example A/C
