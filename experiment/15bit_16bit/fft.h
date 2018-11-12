#pragma once

#include <ostream>
#include <vector>
#include "base.h"

class Fft {
 public:
  enum class Radix {
    Two,
    Three,
    Five,
  };
  enum class Direction {
    Forward,
    Backward,
  };

  void run(const Direction, std::vector<double>& data) const;

 protected:
  Fft(const int64 n,
      const int64 log2n,
      const Radix radix,
      std::vector<double>& work,
      std::vector<double>& table);

 private:
  void Radix2(const int width,
              const int height,
              double* ptr,
              double* x,
              double* y) const;
  void Radix3(const int width,
              const int height,
              double* ptr,
              double* x,
              double* y) const;
  void Radix5(const int width,
              const int height,
              double* ptr,
              double* x,
              double* y) const;

  const int64 n;
  const int64 log2n;
  const Radix radix;
  std::vector<double>& work_;
  std::vector<double>& table_;
};

class Rft : public Fft {
 public:
  Rft(const int64 n,
      const int64 log2n,
      const Radix radix,
      std::vector<double>& work,
      std::vector<double>& table)
      : Fft(n / 2, log2n - 1, radix, work, table), n(n) {}
  void run(const Direction, std::vector<double>& data) const;

 private:
  const int64 n;
};

std::ostream& operator<<(std::ostream&, const Fft::Radix& r);
