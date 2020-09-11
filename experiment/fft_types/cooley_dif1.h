#pragma once

#include <cmath>
#include <vector>

#include "fft.h"

class CooleyDIF1 final : public FFT {
 public:
  CooleyDIF1() = default;
  ~CooleyDIF1() override = default;
  const char* name() const override { return "CooleyDIF1"; }

  void setUp(int n_) override {
    FFT::setUp(n_);

    irevs.resize(n);
    irevs[0] = 0;
    irevs[n - 1] = n - 1;
    int rev = 0;
    for (int i = 1; i < n - 1; ++i) {
      for (int b = n >> 1; b > (rev ^= b); b >>= 1)
        ;
      irevs[i] = rev;
    }

    const double theta = -2 * M_PI / n;
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int k = 0; k < m; ++k) {
        const double t = theta * l * k;
        Complex w{std::cos(t), std::sin(t)};
        ws.push_back(w);
      }
    }
  };

  void tearDown() override { ws.clear(); }

  void dft(Complex* a) override {
    auto pw = ws.data();
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int j = 0; j < l; ++j) {
        for (int k = 0; k < m; ++k) {
          Complex w = pw[k];
          int k0 = 2 * j * m + k;
          int k1 = k0 + m;
          Complex a0 = a[k0];
          Complex a1 = a[k1];
          a[k0] = a0 + a1;
          a[k1] = (a0 - a1) * w;
        }
      }
      pw += m;
    }
  }

  void idft(Complex* a) override {
    auto pw = ws.data() + ws.size();
    for (int l = n / 2, m = 1; l >= 1; l /= 2, m *= 2) {
      pw -= m;
      for (int j = 0; j < l; ++j) {
        for (int k = 0; k < m; ++k) {
          Complex w = pw[k].conj();
          int k0 = 2 * j * m + k;
          int k1 = k0 + m;
          Complex a0 = a[k0];
          Complex a1 = a[k1] * w;
          a[k0] = a0 + a1;
          a[k1] = a0 - a1;
        }
      }
    }
    double inv = 1.0 / n;
    for (int i = 0; i < n; ++i)
      a[i] *= inv;
  }

  void rft(double* x) override {
    Complex* a = reinterpret_cast<Complex*>(x);
    dft(a);

    double x0 = x[0];
    double x1 = x[1];
    x[0] = x0 + x1;
    x[1] = x0 - x1;

    auto f = [](Complex& ai, Complex& aj, const Complex& w) {
      aj = aj.conj();
      Complex c{1 - w.imag, w.real};
      Complex z = c * (ai - aj) * 0.5;
      ai = ai - z;
      aj = (aj + z).conj();
    };

    double th = -M_PI / n;
    Complex w{std::cos(th), std::sin(th)};
    f(a[1], a[n - 1], w);
    for (int i = 2, j = n - 2; i < j; i += 2, j -= 2) {
      Complex w0 = ws[i / 2];
      // TODO: j=>irevs[j]
      f(a[i], a[j], w0);
      f(a[i + 1], a[j - 1], w0 * w);
    }
    a[n / 2] = a[n / 2].conj();
  }

  void irft(double* x) override {
    Complex* a = reinterpret_cast<Complex*>(x);

    double x0 = x[0], x1 = x[1];
    x[0] = (x0 + x1) * 0.5;
    x[1] = (x0 - x1) * 0.5;

    auto f = [](Complex& ai, Complex& aj, const Complex& w) {
      aj = aj.conj();
      Complex c{1 + w.imag, -w.real};
      Complex z = c * (ai - aj) * 0.5;
      ai = ai - z;
      aj = (aj + z).conj();
    };

    double th = M_PI / n;
    Complex w{std::cos(th), std::sin(th)};
    f(a[1], a[n - 1], w);
    for (int i = 2, j = n - 2; i < j; i += 2, j -= 2) {
      Complex w0 = ws[i / 2].conj();
      // TODO: j=>irevs[j]
      f(a[i], a[j], w0);
      f(a[i + 1], a[j - 1], w0 * w);
    }
    a[n / 2] = a[n / 2].conj();
    idft(a);
  }

  std::vector<int> irevs;
  std::vector<Complex> ws;
};
