#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "complex.h"

#define LOG std::cerr << "(" << __LINE__ << "): "

class FMT {
 public:
  FMT(int logn)
      : logn_(logn), log2n_(logn % 2), log4n_(logn / 2), n_(1 << logn) {}
  virtual ~FMT() = default;

  void Convolute(Complex*, Complex*) const;
  bool Test() const;

  // Returns the number of real number operations
  double GetFlops() const { return 8.5 * n_ * log4n_ + 5.0 * n_ * log2n_; }
  virtual double GetMemory() const = 0;

 protected:
  virtual void Rft(Complex*) const = 0;
  virtual void IRft(Complex*) const = 0;

  const int logn_;
  const int log2n_;
  const int log4n_;
  const int n_;
};

class FMTFactoryBase {
 public:
  virtual ~FMTFactoryBase() = default;
  double MeasureConvPerf(int logn) const;

  virtual std::unique_ptr<const FMT> Create(int logn) const = 0;
  virtual const char* name() const = 0;
};

template <typename T>
class FMTFactory : public FMTFactoryBase {
 public:
  ~FMTFactory<T>() override = default;

  static std::unique_ptr<FMTFactoryBase> GetFactory() {
    return std::unique_ptr<FMTFactoryBase>(new FMTFactory<T>());
  }

  std::unique_ptr<const FMT> Create(int logn) const override {
    return std::unique_ptr<const FMT>(new T(logn));
  }

  const char* name() const override { return T::name(); }
};
