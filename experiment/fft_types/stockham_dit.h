#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class StockhamDIT final : public FFT {
 public:
  StockhamDIT(int64 log2k)
    : FFT(log2k) {
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
    for (int64 l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      dft2<backward>(x, y, l, m, pw);
      pw += m;
      std::swap(x, y);
    }
    if (x != a) {
      for (int64 i = 0; i < n; ++i) {
        a[i] = x[i];
      }
    }
    if (backward) {
      double inv = 1.0 / n;
      for (int i = 0; i < n; ++i) {
        a[i] *= inv;
      }
    }
  };

  template <bool backward>
  void dft2(Complex* x, Complex* y, const int64 l, const int64 m, const Complex* pw) const {
    for (int64 k = 0; k < m; ++k) {
      Complex w = backward ? pw[k].conj() : pw[k];
      for (int64 j = 0; j < l; ++j) {
        int64 ix0 = k * l + j;
        int64 ix1 = ix0 + l * m;
        int64 iy0 = k * 2 * l + j;
        int64 iy1 = iy0 + l;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1];
        y[iy0] = x0 + x1;
        y[iy1] = (x0 - x1) * w;
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
  for (int64 i = 0; i < log2n; ++i) {
    m /= 2;
    for (int k = 0; k < m; ++k) {
      const double t = theta * l * k;
      ws.push_back(Complex {std::cos(t), std::sin(t)});
    }
    l *= 2;
  }
}
