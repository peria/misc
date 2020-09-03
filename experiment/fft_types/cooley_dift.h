#pragma once

#include <cmath>
#include <vector>

#include "fft.h"

class CooleyDIFT final : public FFT {
 public:
  CooleyDIFT() = default;
  ~CooleyDIFT() override = default;
  const char* name() const override { return "CooleyDIFT"; }

  void setUp(int n_) override {
    n = n_;
  };

  void dft(Complex* a) override {
    // DIF
    const double theta = -2 * M_PI / n;
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int k = 0; k < m; ++k) {
	const double t = theta * l * k;
	Complex w {std::cos(t), std::sin(t)};
	for (int j = 0; j < l; ++j) {
	  int k0 = 2 * j * m + k;
	  int k1 = k0 + m;
	  Complex a0 = a[k0];
	  Complex a1 = a[k1];
	  a[k0] = a0 + a1;
	  a[k1] = (a0 - a1) * w;
	}
      }
    }
  }

  void idft(Complex* a) override {
    // DIT
    const double theta = 2 * M_PI / n;
    for (int l = n / 2, m = 1; l >= 1; l /= 2, m *= 2) {
      for (int k = 0; k < m; ++k) {
	const double t = theta * l * k;
	Complex w {std::cos(t), std::sin(t)};
	for (int j = 0; j < l; ++j) {
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

  int n;
};
