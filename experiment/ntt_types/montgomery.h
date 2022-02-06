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

template <uint64 P>
struct ReduceHelper {
  template <uint64 bit>
  static inline uint64 reduce(const uint128& x) {
    uint64 y = (x - uint128(uint64(x) * kIP) * P) >> 64;
    return (int64(y) < 0) ? (y + P) : y;
  }
  template <uint64 bit>
  static inline uint64 reduce(uint64 x) {
    if (x == 0) {
      return 0;
    }
    uint64 y = (x - uint128(x * kIP) * P) >> 64;
    return (int64(y) < 0) ? (y + P) : y;
  }

  template <>
  static inline uint64 reduce<1ULL>(const uint128& x) {
    uint128 y = uint128(uint64(x) * kIP) * P;
    uint64 z = (x - y) >> 64;
    return (x < y) ? z + P : z;
  }
  template <>
  static inline uint64 reduce<1ULL>(uint64 x) {
    if (x == 0) {
      return 0;
    }
    uint128 y = uint128(x * kIP) * P;
    uint64 z = (x - y) >> 64;
    return (x < y) ? z + P : z;
  }

 private:
  static constexpr uint64 kIP = ModTrait<P>::kIP;
};

template <uint64 P>
uint64 reduce(const uint128& x) {
  return ReduceHelper<P>::template reduce<(P >> 63)>(x);
}

template <uint64 P>
uint64 reduce(uint64 x) {
  return ReduceHelper<P>::template reduce<(P >> 63)>(x);
}

}  // namespace mod_util

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

  operator uint64() const { return mod_util::reduce<P>(x_); }

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
    x_ = mod_util::reduce<P>(uint128(x_) * y.x_);
    return *this;
  }
  Mod<P> operator*(const Mod<P>& y) {
    Mod<P> p;
    p.x_ = mod_util::reduce<P>(uint128(x_) * y.x_);
    return p;
  }

 private:
  static uint64 init(uint64 x) {
    static constexpr uint64 kR2 = ModTrait<P>::kR2;
    return mod_util::reduce<P>(uint128(kR2) * x);
  }

  static uint64 sub(uint64 a, uint64 b) {
    uint64 c = a - b;
    if (a < b) {
      c += P;
    }
    return c;
  }

  static uint64 add(uint64 a, uint64 b) { return sub(a, P - b); }

 public:
  uint64 x_;
};
