#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class StockhamDIF final : public FFT {
 public:
  StockhamDIF(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~StockhamDIF() override = default;

  void dft(Complex* x, bool backward) const override {
    if (backward)
      dft<true>(x);
    else
      dft<false>(x);
  }

  static const char* name() { return "StockDIF"; }

 private:
  void init();

  template <bool backward>
  void dft(Complex* a) const {
    Complex* x = a;
    Complex* y = const_cast<Complex*>(work.data());
    const Complex* pw = ws.data();
    int64 l = n_;
    int64 m = 1;
    if (log4n_ % 2) {
      std::swap(x, y);
    }
    if (log2n_) {
      l /= 2;
      dft2(a, x, l);
      m *= 2;
    }
    {
      l /= 4;
      dft4<backward>(log2n_ ? x : a, y, l, m, pw);
      pw += m * 3;
      m *= 4;
      std::swap(x, y);
    }
    for (int64 i = 1; i < log4n_; ++i) {
      l /= 4;
      dft4<backward>(x, y, l, m, pw);
      pw += m * 3;
      m *= 4;
      std::swap(x, y);
    }
    if (backward) {
      double inv = 1.0 / n_;
      for (int i = 0; i < n_; ++i) {
        a[i] *= inv;
      }
    }
  };

  void dft2(Complex* x, Complex* y, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 i0 = j;
      int64 i1 = j + l;
      Complex x0 = x[i0];
      Complex x1 = x[i1];
      y[i0] = x0 + x1;
      y[i1] = x0 - x1;
    }
  }

  template <bool backward>
  void dft4(Complex* x, Complex* y, const int64 l, const int64 m, const Complex* pw) const {
    {
      for (int64 j = 0; j < l; ++j) {
        int64 ix0 = j;
        int64 ix1 = l + j;
        int64 ix2 = l * 2 + j;
        int64 ix3 = l * 3 + j;
        int64 iy0 = j;
        int64 iy1 = l * m + j;
        int64 iy2 = l * m * 2 + j;
        int64 iy3 = l * m * 3 + j;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1];
        Complex x2 = x[ix2];
        Complex x3 = x[ix3];
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x1 - x3).i();
        y[iy0] = b0 + b1;
        y[iy1] = backward ? (b2 + b3) : (b2 - b3);
        y[iy2] = b0 - b1;
        y[iy3] = backward ? (b2 - b3) : (b2 + b3);
      }
    }
    for (int64 k = 1; k < m; ++k) {
      Complex w1 = backward ? pw[3 * k].conj() : pw[3 * k];
      Complex w2 = backward ? pw[3 * k + 1].conj() : pw[3 * k + 1];
      Complex w3 = backward ? pw[3 * k + 2].conj() : pw[3 * k + 2];
      for (int64 j = 0; j < l; ++j) {
        int64 ix0 = l * 4 * k + j;
        int64 ix1 = l * (4 * k + 1) + j;
        int64 ix2 = l * (4 * k + 2) + j;
        int64 ix3 = l * (4 * k + 3) + j;
        int64 iy0 = l * k + j;
        int64 iy1 = l * (k + m) + j;
        int64 iy2 = l * (k + m * 2) + j;
        int64 iy3 = l * (k + m * 3) + j;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1] * w1;
        Complex x2 = x[ix2] * w2;
        Complex x3 = x[ix3] * w3;
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x1 - x3).i();
        y[iy0] = b0 + b1;
        y[iy1] = backward ? (b2 + b3) : (b2 - b3);
        y[iy2] = b0 - b1;
        y[iy3] = backward ? (b2 - b3) : (b2 + b3);
      }
    }
  }

  std::vector<Complex> ws;
  std::vector<Complex> work;
};

void StockhamDIF::init() {
  work.resize(n_);
  int64 l = n_;
  int64 m = 1;
  const double theta = -2 * M_PI / n_;
  if (log2n_) {
    l /= 2;
    m *= 2;
  }
  for (int64 i = 0; i < log4n_; ++i) {
    l /= 4;
    for (int64 k = 0; k < m; ++k) {
      double t = theta * l * k;
      ws.push_back(Complex {std::cos(t), std::sin(t)});
      ws.push_back(Complex {std::cos(2*t), std::sin(2*t)});
      ws.push_back(Complex {std::cos(3*t), std::sin(3*t)});
    }
    m *= 4;
  }
};
