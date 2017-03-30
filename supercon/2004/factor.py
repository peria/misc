#!/usr/bin/python

import fractions
import time
from itertools import chain, tee, dropwhile

# Prime generator is from
# http://qiita.com/ttatsf/items/da3bc1a26740d63586b4#_reference-2ac8e0cb2fe5249a2c93

def primes():
  cropped_primables = (lambda s, p :
      [x + y for x in xrange(s, p * p - 5, 6) for y in [0, 4]]
  )

  zoned_primes = (lambda m, s, p :
      [n for n in cropped_primables(s, p) if fractions.gcd(m, n) == 1]
  )

  primes_list = [2, 3, 5]
  last_prime = primes_list[-1]
  primes_iter = iter(primes_list)

  seeds_iter = iter([])
  m, s, p = 1, 7, 5

  while True:
    prime = primes_iter.next()
    yield prime

    if prime == last_prime :
      primes_list = zoned_primes(m, s, p)
      last_prime = primes_list[-1]

      it1, it2 = tee(primes_list)
      primes_iter = it1
      seeds_iter = chain(seeds_iter, it2)

      m = m * p
      s = p * p
      p = seeds_iter.next()


def factor(n):
  due_time = time.time() + 180  # 180 secs
  # Use p-1 method
  k, a = 0, 2
  ps = primes()
  while True:
    p = ps.next()
    if p < 10000:
      for _ in xrange(5):
        a = pow(a, p, n)
      continue
    a = pow(a, p, n)
    k = k + 1
    if k % 100000 != 0:
      continue
    g = fractions.gcd(a - 1, n)
    if 1 < g:
      return g
    if due_time < time.time():
      break
  return fractions.gcd(a - 1, n)


def test(n):
  start = time.time()
  p = factor(n)
  end = time.time()

  if p == n:
    print 'Too many primes to factor %x' % n
  elif p == 1:
    print 'Fail to factor %x in 180 sec.' % n
  else:
    print '%x = %x * %x' % (n, p, n / p)
    print '%d = %d * %d' % (n, p, n / p)
  print 'Time: %f sec' % (end - start)


if __name__ == '__main__':
  # These numbers are in http://www.gsic.titech.ac.jp/supercon/supercon2004/kaisetsu/
  test(0xd44cdd3dc3436e0c04cdedb0f)  # Contest
  test(0x851bed3b13e22bd06bc3a0ba7)  # Example A/C
