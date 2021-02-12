#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class PMP : public FFT {
 public:
  PMP(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~PMP() override = default;

  static const char* name() { return "PMP"; }

  void dft(Complex* x, bool backward) const override {
    if (backward)
      idft(x);
    else
      dft(x);
  }

 private:
  void init();

  void dft(Complex* a) const {
    auto* pw = ws.data();
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log4n; ++i) {
      m /= 4;
      dft4(a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 i = 0; i < log2n; ++i) {
      m /= 2;
      dft2(a, l);
      l *= 2;
    }
  }

  void idft(Complex* a) const {
    auto* pw = ws.data() + ws.size();
    int64 l = n;
    int64 m = 1;
    if (log2n) {
      l /= 2;
      dft2(a, l);
      m *= 2;
    }
    for (int64 i = 0; i < log4n; ++i) {
      l /= 4;
      pw -= m * 3;
      idft4(a, l, m, pw);
      m *= 4;
    }
    double inv = 1.0 / n;
    for (int64 i = 0; i < n; ++i)
      a[i] *= inv;
  }

  void dft2(Complex* a, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 k0 = 2 * j;
      int64 k1 = 2 * j + 1;
      Complex a0 = a[k0];
      Complex a1 = a[k1];
      a[k0] = a0 + a1;
      a[k1] = a0 - a1;
    }
  }

  void dft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 k0 = 4 * j * m;
        int64 k1 = 4 * j * m + m;
        int64 k2 = 4 * j * m + 2 * m;
        int64 k3 = 4 * j * m + 3 * m;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[k0] = b0 + b1;
        a[k1] = b0 - b1;
        a[k2] = b2 + b3;
        a[k3] = b2 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k * 3];
        Complex w2 = ws[k * 3 + 1];
        Complex w3 = ws[k * 3 + 2];
        int64 k0 = 4 * j * m + k;
        int64 k1 = 4 * j * m + m + k;
        int64 k2 = 4 * j * m + 2 * m + k;
        int64 k3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[k0] = b0 + b1;
        a[k1] = (b0 - b1) * w2;
        a[k2] = (b2 + b3) * w1;
        a[k3] = (b2 - b3) * w3;
      }
    }
  }

  void idft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 k0 = 4 * j * m;
        int64 k1 = 4 * j * m + m;
        int64 k2 = 4 * j * m + 2 * m;
        int64 k3 = 4 * j * m + 3 * m;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[k0] = b0 + b2;
        a[k1] = b1 + b3;
        a[k2] = b0 - b2;
        a[k3] = b1 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[3 * k].conj();
        Complex w2 = ws[3 * k + 1].conj();
        Complex w3 = ws[3 * k + 2].conj();
        int64 k0 = 4 * j * m + k;
        int64 k1 = 4 * j * m + m + k;
        int64 k2 = 4 * j * m + 2 * m + k;
        int64 k3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[k0];
        Complex a1 = a[k1] * w2;
        Complex a2 = a[k2] * w1;
        Complex a3 = a[k3] * w3;
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[k0] = b0 + b2;
        a[k1] = b1 + b3;
        a[k2] = b0 - b2;
        a[k3] = b1 - b3;
      }
    }
  }

  std::vector<Complex> ws;
};

void PMP::init() {
  int64 l = 1;
  int64 m = n;
  const double theta = -2 * M_PI / n;
  for (int64 i = 0; i < log4n; ++i) {
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
