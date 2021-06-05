#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "complex.h"
#include "fft.h"

// PMP with DFS algorthm. Focus on in-cache area.
class DPMP : public FFT {
 public:
  DPMP(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~DPMP() override = default;

  static const char* name() { return "DPMP"; }

  void dft(Complex* x, bool backward) const override {
    if (backward)
      idft(x);
    else
      dft(x);
  }

  void rft(Complex* x, bool backward) const override {
    if (backward) {
      dft(x, backward);
      const Complex& wq1 = wq.conj();
      const Complex& wq2 = wq1 * wq1;
      const Complex& wq3 = wq2 * wq1;
      for (int i = 0; i < n_ / 4; ++i) {
        const Complex& w = ws[3 * i].conj();
        x[4 * i] *= w;
        x[4 * i + 1] *= w * wq1;
        x[4 * i + 2] *= w * wq2;
        x[4 * i + 3] *= w * wq3;
      }
    } else {
      const Complex& wq1 = wq;
      const Complex& wq2 = wq1 * wq1;
      const Complex& wq3 = wq2 * wq1;
      for (int i = 0; i < n_ / 4; ++i) {
        const Complex& w = ws[3 * i];
        x[4 * i] *= w;
        x[4 * i + 1] *= w * wq1;
        x[4 * i + 2] *= w * wq2;
        x[4 * i + 3] *= w * wq3;
      }
      dft(x, backward);
    }
  }

 private:
  void init();

  void dft(Complex* a) const {
    auto* pw = ws.data();
    int64 l = 1;
    int64 m = n_;
    int64 log4n_local = std::min<int64>(log4n_, 8);
    int64 log4n_all = log4n_ - log4n_local;
    for (int64 i = 0; i < log4n_all; ++i) {
      m /= 4;
      dft4(a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 k = 0; k < l; ++k) {
      Complex* al = a + k * m;
      auto* pwl = pw;
      int64 ll = 1;
      int64 ml = m;
      for (int64 i = 0; i < log4n_local; ++i) {
        ml /= 4;
        dft4(al, ll, ml, pwl);
        pwl += ml * 3;
        ll *= 4;
      }
      for (int64 i = 0; i < log2n_; ++i) {
        ml /= 2;
        dft2(al, ll);
        ll *= 2;
      }
    }
  }

  void idft(Complex* a) const {
    const int64 log4n_local = std::min<int64>(log4n_, 8);
    const int64 log4n_all = log4n_ - log4n_local;
    const int64 ll = 1LL << (log2n_ + log4n_local * 2);
    const int64 ml = 1LL << (log4n_all * 2);
    for (int64 k = 0; k < ml; ++k) {
      auto* pw = ws.data() + ws.size();
      Complex* al = a + k * ll;
      int64 l = ll;
      int64 m = 1;
      if (log2n_) {
        l /= 2;
        dft2(a, l);
        m *= 2;
      }
      for (int64 i = 0; i < log4n_; ++i) {
        l /= 4;
        pw -= m * 3;
        idft4(a, l, m, pw);
        m *= 4;
      }
    }

    if (log4n_all) {
      auto* pw = ws.data() + ws.size();
      int64 l = n_;
      int64 m = 1;
      if (log2n_) {
        l /= 2;
        m *= 2;
      }
      for (int64 i = 0; i < log4n_; ++i) {
        l /= 4;
        pw -= m * 3;
        m *= 4;
      }
      for (int64 i = 0; i < log4n_all; ++i) {
        l /= 4;
        pw -= m * 3;
        idft4(a, l, m, pw);
        m *= 4;
      }
    }
    double inv = 1.0 / n_;
    for (int64 i = 0; i < n_; ++i)
      a[i] *= inv;
  }

  void dft2(Complex* a, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 i0 = 2 * j;
      int64 i1 = 2 * j + 1;
      Complex a0 = a[i0];
      Complex a1 = a[i1];
      a[i0] = a0 + a1;
      a[i1] = a0 - a1;
    }
  }

  void dft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[i0] = b0 + b1;
        a[i1] = b0 - b1;
        a[i2] = b2 + b3;
        a[i3] = b2 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k * 3];
        Complex w2 = ws[k * 3 + 1];
        Complex w3 = ws[k * 3 + 2];
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[i0] = b0 + b1;
        a[i1] = (b0 - b1) * w2;
        a[i2] = (b2 + b3) * w1;
        a[i3] = (b2 - b3) * w3;
      }
    }
  }

  void idft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[i0] = b0 + b2;
        a[i1] = b1 + b3;
        a[i2] = b0 - b2;
        a[i3] = b1 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[3 * k].conj();
        Complex w2 = ws[3 * k + 1].conj();
        Complex w3 = ws[3 * k + 2].conj();
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[i0];
        Complex a1 = a[i1] * w2;
        Complex a2 = a[i2] * w1;
        Complex a3 = a[i3] * w3;
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[i0] = b0 + b2;
        a[i1] = b1 + b3;
        a[i2] = b0 - b2;
        a[i3] = b1 - b3;
      }
    }
  }

  std::vector<Complex> ws;
  Complex wq;
};

void DPMP::init() {
  int64 l = 1;
  int64 m = n_;
  const double theta = -2 * M_PI / n_;
  wq = Complex {std::cos(theta / 4), -std::sin(theta / 4)};
  for (int64 i = 0; i < log4n_; ++i) {
    m /= 4;
    for (int64 k = 0; k < m; ++k) {
      double t = theta * l * k;
      ws.push_back(Complex {std::cos(t), std::sin(t)});
      ws.push_back(Complex {std::cos(2*t), std::sin(2*t)});
      ws.push_back(Complex {std::cos(3*t), std::sin(3*t)});
    }
    l *= 4;
  }
}
