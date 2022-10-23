#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "complex.h"

#define LOG std::cerr << "(" << __LINE__ << "): "

class FMT {
 public:
  enum class Type { kFFT, kNTT };

  FMT(int log2n, int log4n)
      : log2n_(log2n), log4n_(log4n), n_(1 << (log2n_ + log4n_ * 2)) {}
  virtual ~FMT() = default;
  virtual Type type() const = 0;

  template <typename T>
  const T& As() const {
    return *dynamic_cast<const T*>(this);
  }

 protected:
  const int log2n_;
  const int log4n_;
  const int n_;
};

class FFT : public FMT {
 public:
  FFT(int log2n, int log4n = 0) : FMT(log2n, log4n) {}
  Type type() const override { return Type::kFFT; }

  void Convolution(std::vector<Complex>& x, std::vector<Complex>& y) const;

 protected:
  virtual void Dft(Complex*) const = 0;
  virtual void IDft(Complex*) const = 0;
  virtual void Rft(Complex*) const = 0;
  virtual void IRft(Complex*) const = 0;
};

class NTT : public FMT {
  Type type() const override { return Type::kNTT; }
};

class FMTFactory {
 public:
  virtual ~FMTFactory() = default;
  // logn is log of length of integer
  virtual std::shared_ptr<FMT> Create(int logn) const = 0;
  virtual const char* name() const = 0;
};
