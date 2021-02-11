#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class StockhamDIT final : public FFT {
 public:
  StockhamDIT(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~StockhamDIT() override = default;

  void dft(Complex* x, bool backward) const override {
    if (backward)
      dft<true>(x);
    else
      dft<false>(x);
  }

  static const char* name() { return "StockDIT"; }

 private:
  void init();

  template <bool backward>
  void dft(Complex* a) const {
    Complex* x = a;
    Complex* y = const_cast<Complex*>(work.data());
    const Complex* pw = ws.data();
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log4n - 1; ++i) {
      m /= 4;
      dft4<backward>(x, y, l, m, pw);
      pw += m * 3;
      l *= 4;
      std::swap(x, y);
    }
    {
      m /= 4;
      dft4<backward>(x, a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    if (log2n) {
      m /= 2;
      dft2(a, a, l, m);
      pw += m;
      l *= 2;
    }
    if (backward) {
      double inv = 1.0 / n;
      for (int i = 0; i < n; ++i) {
        a[i] *= inv;
      }
    }
  };

  void dft2(Complex* x, Complex* y, const int64 l, const int64 m) const {
    for (int64 k = 0; k < m; ++k) {
      for (int64 j = 0; j < l; ++j) {
        int64 i0 = j;
        int64 i1 = j + l;
        Complex x0 = x[i0];
        Complex x1 = x[i1];
        y[i0] = x0 + x1;
        y[i1] = x0 - x1;
      }
    }
  };

  template <bool backward>
  void dft4(Complex* x, Complex* y, const int64 l, const int64 m, const Complex* pw) const {
    {
      for (int64 j = 0; j < l; ++j) {
        int64 i0 = j;
        int64 i1 = l + j;
        int64 i2 = 2 * l + j;
        int64 i3 = 3 * l + j;
        Complex x0 = x[i0];
        Complex x1 = x[i1];
        Complex x2 = x[i2];
        Complex x3 = x[i3];
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x3 - x1).i();
        y[i0] = b0 + b1;
        y[i1] = b0 - b1;
        y[i2] = b2 + b3;
        y[i3] = b2 - b3;
      }
    }
    for (int64 k = 1; k < m; ++k) {
      Complex w1 = backward ? pw[k * 3].conj() : pw[k * 3];
      Complex w2 = backward ? pw[k * 3 + 1].conj() : pw[k * 3 + 1];
      Complex w3 = backward ? pw[k * 3 + 2].conj() : pw[k * 3 + 2];
      for (int64 j = 0; j < l; ++j) {
        int64 ix0 = k * l + j;
        int64 ix1 = (k + m) * l + j;
        int64 ix2 = (k + m * 2) * l + j;
        int64 ix3 = (k + m * 3) * l + j;
        int64 iy0 = k * 4 * l + j;
        int64 iy1 = (k * 4 + 1) * l + j;
        int64 iy2 = (k * 4 + 2) * l + j;
        int64 iy3 = (k * 4 + 3) * l + j;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1];
        Complex x2 = x[ix2];
        Complex x3 = x[ix3];
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x3 - x1).i();
        y[iy0] = b0 + b1;
        y[iy1] = (b0 - b1) * w2;
        y[iy2] = (b2 + b3) * w1;
        y[iy3] = (b2 - b3) * w3;
      }
    }
  };

  std::vector<Complex> ws;
  std::vector<Complex> work;
};

void StockhamDIT::init() {
  work.resize(n);
  int64 l = 1;
  int64 m = n;
  const double theta = -2 * M_PI / n;
  for (int64 i = 0; i < log4n; ++i) {
    m /= 4;
    for (int k = 0; k < m; ++k) {
      const double t = theta * l * k;
      ws.push_back(Complex {std::cos(t), std::sin(t)});
      ws.push_back(Complex {std::cos(2*t), std::sin(2*t)});
      ws.push_back(Complex {std::cos(3*t), std::sin(3*t)});
    }
    l *= 4;
  }
}
