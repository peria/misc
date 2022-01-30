#pragma once

#include "base.h"

namespace mod_util {
uint64 inverse(uint64 x) {
  uint64 inv = x;
  for (int64 i = 0; i < 5; ++i) {
    inv *= 2 - inv * x;
  }
  return inv;
}
uint64 neg_inverse(uint64 x) {
  uint64 inv = 0;
  uint64 s = 0;
  for (uint64 b = 1; b; b <<= 1) {
    if (s & b)
      continue;
    s += b * x;
    inv += b;
  }
  return inv;
}
}  // namespace mod_util

template <uint64 P>
struct ModTrait {
  static const uint64 kIP = mod_util::inverse(P);
  static constexpr uint64 kR2 = uint128(-P) % P;
};

template <>
struct ModTrait<0xffffffff00000001ULL> {
  static constexpr uint64 kIP = 0x100000001ULL;
  static constexpr uint64 kR2 = 0xfffffffe00000001ULL;
};

template <uint64 P>
class Mod {
 public:
  Mod() = default;
  explicit Mod(uint64 y) : x_(init(y)) {}
  Mod<P>& operator=(const Mod<P>& y) {
    x_ = y.x_;
    return *this;
  }
  Mod<P>& operator=(uint64 y) {
    x_ = init(y);
    return *this;
  }

  operator uint64() const { return reduce(x_); }

  Mod<P> operator+(const Mod<P>& y) {
    Mod<P> s;
    s.x_ = add(x_, y.x_);
    return s;
  }
  Mod<P>& operator+=(const Mod<P>& y) {
    x_ = add(x_, y.x_);
    return *this;
  }

  Mod<P> operator-(const Mod<P>& y) {
    Mod<P> s;
    s.x_ = sub(x_, y.x_);
    return s;
  }
  Mod<P>& operator-=(const Mod<P>& y) {
    x_ = sub(x_, y.x_);
    return *this;
  }

  Mod<P>& operator*=(const Mod<P>& y) {
    x_ = reduce(uint128(x_) * y.x_);
    return *this;
  }
  Mod<P> operator*(const Mod<P>& y) {
    Mod<P> p;
    p.x_ = reduce(uint128(x_) * y.x_);
    return p;
  }

 private:
  static uint64 init(uint64 x) {
    static constexpr uint64 kR2 = ModTrait<P>::kR2;
    return reduce(uint128(kR2) * x);
  }

  static uint64 sub(uint64 a, uint64 b) {
    uint64 c = a - b;
    if (a < b) {
      c += P;
    }
    return c;
  }

  static uint64 add(uint64 a, uint64 b) { return sub(a, P - b); }

  static uint64 reduce(uint64 x) {
    static constexpr uint64 kIP = ModTrait<P>::kIP;
    if (x == 0) {
      return 0;
    }
    uint128 y = uint128(x * kIP) * P;
    uint64 z = (x - y) >> 64;
    return (x < y) ? z + P : z;
  }

  static uint64 reduce(const uint128& x) {
    static constexpr uint64 kIP = ModTrait<P>::kIP;
    uint128 y = uint128(uint64(x) * kIP) * P;
    uint64 z = (x - y) >> 64;
    return (x < y) ? z + P : z;
  }

 public:
  uint64 x_;
};
