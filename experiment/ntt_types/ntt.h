#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include "base.h"
#include "mod.h"

class NTT {
 public:
  static constexpr uint64 P = 0xffffffff00000001ULL;
  static constexpr uint64 W = 7ULL;
  static constexpr uint64 IW = 0x249249246db6db6eULL;
  using ElementType = Mod<P>;

  NTT(int64 log2n, int64 log4n = 0)
      : log2n_(log2n),
        log4n_(log4n),
        logn_(log2n_ + log4n_ * 2),
        n_(1LL << (logn_)) {}

  virtual ~NTT() = default;
  virtual void ntt(ElementType* x, bool backward) const = 0;

  double getIop() const { return 1.5 * n_ * logn_; }

  int64 size() const { return n_; }

 protected:
  const int64 log2n_ = 0;
  const int64 log4n_ = 0;
  const int64 logn_ = 0;
  const int64 n_ = 1;
};

class NTTFactoryBase {
 public:
  NTTFactoryBase() = default;
  virtual ~NTTFactoryBase() = default;

  virtual NTT* Create(int64 logn) = 0;
  virtual const char* name() const = 0;
};

template <typename T>
class NTTFactory final : public NTTFactoryBase {
 public:
  NTTFactory() = default;
  ~NTTFactory() override = default;

  NTT* Create(int64 logn) override { return new T(logn); }
  const char* name() const override { return T::name(); }
};
