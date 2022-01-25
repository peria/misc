#pragma once

#include "base.h"
#include "mod_base.h"

template <>
class Mod<0xffffffff00000001ULL> {
  static constexpr uint64 P = 0xffffffff00000001ULL;

 public:
  Mod() = default;
  explicit Mod(uint64 x) : x_(x) {}
  Mod<P>& operator=(const Mod<P>& x) {
    x_ = x;
    return *this;
  }
  Mod<P>& operator=(uint64 x) {
    x_ = x;
    return *this;
  }

  operator uint64() const { return x_; }

  Mod<P> operator+(const Mod<P>& y) { return operator-(Mod<P>(P - y.x_)); }
  Mod<P>& operator+=(const Mod<P>& y) { return operator-=(Mod<P>(P - y.x_)); }

  Mod<P> operator-(const Mod<P>& y) { return Mod<P>(sub(x_, y.x_)); }
  Mod<P>& operator-=(const Mod<P>& y) {
    x_ = sub(x_, y.x_);
    return *this;
  }

  Mod<P>& operator*=(const Mod<P>& y) {
    uint128 prod = uint128(x_) * y.x_;
    uint64 h = prod >> 64;
    uint64 l = prod;
    x_ = reduce(l, h);
    return *this;
  }
  Mod<P> operator*(const Mod<P>& y) {
    uint128 prod = uint128(x_) * y.x_;
    uint64 h = prod >> 64;
    uint64 l = prod;
    return Mod<P>(reduce(l, h));
  }

  Mod<P>& operator>>=(int64 k) {
    x_ >>= k;
    return *this;
  }

 private:
  static uint64 reduce(uint64 a0, uint64 a1) {
    static constexpr uint64 kLowerMask = 0xffffffffULL;
    uint64 a10 = a1 & kLowerMask;
    uint64 a11 = a1 >> 32;
    if (a0 >= P)
      a0 -= P;
    a0 = sub(a0, a10);
    a0 = sub(a0, a11);
    return add(a0, a10 << 32);
  }

  static uint64 sub(uint64 x, uint64 y) {
    uint64 d = x - y;
    if (x < y) {
      d += P;
    }
    return d;
  }

  static uint64 add(uint64 x, uint64 y) { return sub(x, P - y); }

  uint64 x_;
};
