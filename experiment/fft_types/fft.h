#pragma once

#include <cmath>
#include "complex.h"

class FFT {
 public:
  FFT(int64 log2n, int64 log4n=0)
    : log2n(log2n),
      log4n(log4n),
      logn(log2n + log4n * 2),
      n(1LL << (log2n + log4n * 2)) {}
      
  virtual ~FFT() = default;
  virtual void dft(Complex* x, bool backward) const = 0;
  virtual double getFlop() const = 0;

 protected:
  const int64 log2n = 0;
  const int64 log4n = 0;
  const int64 logn = 0;
  const int64 n = 1;
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
