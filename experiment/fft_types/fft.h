#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include "complex.h"

class FFT {
 public:
  FFT(int64 log2n, int64 log4n=0)
    : log2n_(log2n),
      log4n_(log4n),
      logn_(log2n_ + log4n_ * 2),
      n_(1LL << (logn_)) {}
      
  virtual ~FFT() = default;
  virtual void dft(Complex* x, bool backward) const = 0;
  virtual void rft(Complex* x, bool backward) const {
    // Not implemented methods
    assert(false);
  }

  double getFlop() const {
    double flop = 34 * log4n_ * n_ / 4 + 4 * log2n_ * n_ / 2;
    return flop;
  }

  int64 size() const { return n_; }

 protected:
  const int64 log2n_ = 0;
  const int64 log4n_ = 0;
  const int64 logn_ = 0;
  const int64 n_ = 1;
};

class FFTFactoryBase {
 public:
  FFTFactoryBase() = default;
  virtual ~FFTFactoryBase() = default;

  virtual FFT* Create(int64 logn) = 0;
  virtual const char* name() const = 0;
};

template <typename T>
class FFTFactory final : public FFTFactoryBase {
 public:
  FFTFactory() = default;
  ~FFTFactory() override = default;

  FFT* Create(int64 logn) override { return new T(logn); }
  const char* name() const override { return T::name(); }
};
