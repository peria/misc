#!/usr/bin/python3

def main():
    R = 2 ** 63
    N = 2 ** 27
    m = 3
    for k in range(R // N, 0, -1):
        p = k * N + 1
        if not is_prime(p):
            continue
        r = get_root(p)
        ir = pow(r, p - 2, p)
        r2 = 2**128 % p
        inv = p
        for _ in range(5):
            inv = inv * (2 - inv * p) % 2 ** 64
        print(hex(p), hex(r), hex(ir), hex(r2), hex(inv))
        m -= 1
        if m == 0:
            break


def get_root(p):
    for a in range(2, p):
        if pow(a, (p - 1) // 2, p) == p - 1:
            return a


def is_prime(x):
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    if any([x % p == 0 for p in primes]):
        return False
    t, k = (x - 1, 0)
    while t % 2 == 0:
        t //= 2
        k += 1
    for a in (2, 3, 5, 7, 11):
        y = pow(a, t, x)
        if y == 1:
            continue
        z = False
        for _ in range(k):
            if y == x - 1:
                z = True
                break
            y = pow(y, 2, x)
        if not z:
            return False
    return True


if __name__ == '__main__':
    main()
