#pragma once

#include "base.h"

template <uint64 P>
class Mod {
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

  Mod<P> operator-(const Mod<P>& y) {
    uint64 d = x_ - y.x_;
    if (x_ < y.x_) {
      d += P;
    }
    return Mod<P>(d);
  }
  Mod<P>& operator-=(const Mod<P>& y) {
    uint64 d = x_ - y.x_;
    if (x_ < y.x_) {
      d += P;
    }
    x_ = d;
    return *this;
  }

  Mod<P>& operator*=(const Mod<P>& y) {
    uint128 prod = uint128(x_) * y.x_;
    x_ = prod % P;
    return *this;
  }
  Mod<P> operator*(const Mod<P>& y) {
    uint128 prod = uint128(x_) * y.x_;
    return Mod<P>(uint64(prod % P));
  }

 private:
  uint64 x_;
};
