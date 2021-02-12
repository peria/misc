#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class Cooley : public FFT {
 public:
  Cooley(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~Cooley() override = default;

  static const char* name() { return "Cooley"; }

  void dft(Complex* x, bool backward) const override {
    if (backward)
      dft<true>(x);
    else
      dft<false>(x);
  }

 private:
  void init();
  void sort(Complex*) const;

  template <bool backward>
  void dft(Complex* a) const {
    auto* pw = ws.data();
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log4n; ++i) {
      m /= 4;
      dft4<backward>(a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 i = 0; i < log2n; ++i) {
      m /= 2;
      dft2(a, l);
      l *= 2;
    }
    sort(a);
    if (backward) {
      double inv = 1.0 / n;
      for (int64 i = 0; i < n; ++i) {
        a[i] *= inv;
      }
    }
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

  template <bool backward=false>
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
        Complex w1 = backward ? ws[k * 3].conj() : ws[k * 3];
        Complex w2 = backward ? ws[k * 3 + 1].conj() : ws[k * 3 + 1];
        Complex w3 = backward ? ws[k * 3 + 2].conj() : ws[k * 3 + 2];
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

  std::vector<Complex> ws;
  std::vector<int64> ip;
};

void Cooley::init() {
  {
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

  {
    ip.resize(1LL << (logn / 2));
    ip[0] = 0;
    int64 l = 1;
    int64 m = n;
    while (2 * l < m) {
      m = m / 2;
      for (int64 j = 0; j < l; ++j) {
        ip[l + j] = ip[j] + m;
      }
      l = l * 2;
    }
  }  
}

void Cooley::sort(Complex* a) const {
  int64 l = 1;
  int64 m = n;
  if (l == m) {
    for (int64 i = 1; i < l; ++i) {
      for (int64 j = 0; j < i; ++j) {
        int64 ji = j + ip[i];
        int64 ij = i + ip[j];
        std::swap(a[ji], a[ij]);
      }
    }
  } else {
    for (int64 i = 1; i < l; ++i) {
      for (int64 j = 0; j < i; ++j) {
        int64 ji = j + ip[i];
        int64 ij = i + ip[j];
        std::swap(a[ji], a[ij]);
        std::swap(a[ji + l], a[ij + l]);
      }
    }
  }
}
