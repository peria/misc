#pragma once

#include <cmath>
#include <vector>

#include "fft.h"

class CooleyDITF final : public FFT {
 public:
  CooleyDITF() = default;
  ~CooleyDITF() override = default;
  const char* name() const override { return "CooleyDITF"; }

  void setUp(int n_) override {
    n = n_;
  };

  void dft(Complex* a) override {
    // DIT
    const double theta = - 2 * M_PI / n;
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int j = 0; j < l; ++j) {
	const double t = theta * bitrev(j, n);
	Complex w {std::cos(t), std::sin(t)};
	for (int k = 0; k < m; ++k) {
	  int k0 = 2 * j * m + k;
	  int k1 = k0 + m;
	  Complex a0 = a[k0];
	  Complex a1 = a[k1] * w;
	  a[k0] = a0 + a1;
	  a[k1] = a0 - a1;
	}
      }
    }
  };

  void idft(Complex* a) override {
    // DIF
    const double theta = 2 * M_PI / n;
    for (int l = n / 2, m = 1; m < n; l /= 2, m *= 2) {
      for (int j = 0; j < l; ++j) {
	const double t = theta * bitrev(j, n);
	Complex w {std::cos(t), std::sin(t)};
	for (int k = 0; k < m; ++k) {
	  int k0 = 2 * j * m + k;
	  int k1 = k0 + m;
	  Complex a0 = a[k0];
	  Complex a1 = a[k1];
	  a[k0] = a0 + a1;
	  a[k1] = (a0 - a1) * w;
	}
      }
    }
    double inv = 1.0 / n;
    for (int i = 0; i < n; ++i)
      a[i] *= inv;
  };

  int bitrev(int j, int n) {
    int ret = 0;
    for(int b = n >> 1; b ; b >>= 1, j >>= 1)
      if (j & 1)
	ret |= b;
    return ret;
  }

  int n;
};
