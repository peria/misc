#include "integer.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "complex.h"

void Integer::Convolution(const FMT& fmt,
                          const int n,
                          std::vector<Digit>& x,
                          const std::vector<Digit>& y) {
  switch (fmt.type()) {
    case FMT::Type::kFFT: {
      std::vector<Complex> fx(2 * n);
      std::vector<Complex> fy(2 * n);
      Split(fx, x);
      Split(fy, y);
      fmt.As<FFT>().Convolution(fx, fy);
      Merge(x, fx);
      break;
    }
    case FMT::Type::kNTT: {
      break;
    }
  }
}

double Integer::GetFlops(int logn) {
  // Here N is the number of complex numbers in case of FFT.
  int n = 1 << (logn + 1);
  int log2n = (logn + 1) % 2;
  int log4n = (logn + 1) / 2;
  return (log2n * 8 + log4n * 8.5) * n;
}

void Integer::Split(std::vector<Complex>& fx, const std::vector<Digit>& x) {
  static constexpr int64 B1 = 10000;
  static constexpr int64 B2 = B1 * B1;
  static constexpr int64 B3 = B2 * B1;
  for (int i = 0; i < fx.size() / 4; ++i) {
    Digit xi = x[i];
    Digit xj = x[i + fx.size() / 2];
    double xi0 = xi % B1;
    double xi1 = xi / B1 % B1;
    double xi2 = xi / B2 % B1;
    double xi3 = xi / B3;
    double xj0 = xj % B1;
    double xj1 = xj / B1 % B1;
    double xj2 = xj / B2 % B1;
    double xj3 = xj / B3;
    fx[4 * i + 0] = Complex(xi0, xj0);
    fx[4 * i + 1] = Complex(xi1, xj1);
    fx[4 * i + 2] = Complex(xi2, xj2);
    fx[4 * i + 3] = Complex(xi3, xj3);
  }
}

void Integer::Merge(std::vector<Digit>& x, const std::vector<Complex>& fx) {
  static constexpr int64 B1 = 10000;
  static constexpr int64 B2 = B1 * B1;
  static constexpr int64 B3 = B2 * B1;
  static constexpr int64 B4 = B2 * B2;
  Digit cr = 0;
  Digit ci = 0;
  const int n = fx.size() / 4;
  for (int i = 0; i < n; ++i) {
    Digit x0r = std::round(fx[4 * i + 0].real);
    Digit x0i = std::round(-fx[4 * i + 0].imag);
    Digit x1r = std::round(fx[4 * i + 1].real);
    Digit x1i = std::round(-fx[4 * i + 1].imag);
    Digit x2r = std::round(fx[4 * i + 2].real);
    Digit x2i = std::round(-fx[4 * i + 2].imag);
    Digit x3r = std::round(fx[4 * i + 3].real);
    Digit x3i = std::round(-fx[4 * i + 3].imag);
    Digit xr = cr + x0r + x1r % B3 * B1 + x2r % B2 * B2 + x3r % B1 * B3;
    cr = x1r / B1 + x2r / B2 + x3r / B3 + xr / B4;
    x[i] = xr % B4;
    Digit xi = ci + x0i + x1i % B3 * B1 + x2i % B2 * B2 + x3i % B1 * B3;
    ci = x1i / B1 + x2i / B2 + x3i / B3 + xi / B4;
    x[i + n] = ci % B4;
  }
  // TODO: Carry `cr` to x[i+n/4]
}
