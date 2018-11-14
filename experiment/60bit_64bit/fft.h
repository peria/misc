#pragma once

#include "base.h"
#include "complex.h"
#include <ostream>
#include <vector>

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

  void run(const Direction, std::vector<double> &data) const;

protected:
  Fft(const int64 n, const int64 log2n, const Radix radix,
      std::vector<double> &work, std::vector<double> &table);

private:
  void setTable(int64 r, const int64 height, Complex *table);
  void radix2(const int height, Complex *ptr, Complex *x, Complex *y) const;
  void radix4(const int width, const int height, Complex *ptr, Complex *x,
              Complex *y) const;
  void radix8(const int width, const int height, Complex *ptr, Complex *x,
              Complex *y) const;
  void radix3(const int width, const int height, Complex *x, Complex *y) const;
  void radix5(const int width, const int height, Complex *x, Complex *y) const;

  const int64 n;
  const int64 log2n;
  const int64 log4n;
  const int64 log8n;
  const Radix radix;
  std::vector<double> &work_;
  std::vector<double> &table_;
};

class Rft : public Fft {
public:
  Rft(const int64 n, const int64 log2n, const Radix radix,
      std::vector<double> &work, std::vector<double> &table)
      : Fft(n / 2, log2n - 1, radix, work, table), n(n) {}
  void run(const Direction, std::vector<double> &data) const;

private:
  const int64 n;
};

std::ostream &operator<<(std::ostream &, const Fft::Radix &r);
