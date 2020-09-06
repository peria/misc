#pragma once

#include <cmath>
#include <vector>

#include "fft.h"

class CooleyDIT1 final : public FFT {
 public:
  CooleyDIT1() = default;
  ~CooleyDIT1() override = default;
  const char* name() const override { return "CooleyDIT1"; }

  void setUp(int n_) override {
    FFT::setUp(n_);

    irevs.resize(n);
    irevs[0] = 0;
    irevs[n - 1] = n - 1;
    int rev = 0;
    for (int i = 1; i < n - 1; ++i) {
      for (int b = n >> 1; b > (rev ^= b); b >>= 1);
      irevs[i] = rev;
    }

    const double theta = -2 * M_PI / n;
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int j = 0; j < l; ++j) {
	const double t = theta * irevs[j * 2];
	Complex w {std::cos(t), std::sin(t)};
	ws.push_back(w);
      }
    }
  };

  void tearDown() override {
    ws.clear();
  }

  void dft(Complex* a) override {
    auto iw = ws.begin();
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int j = 0; j < l; ++j) {
	Complex w = *iw++;
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
    for (int i = 1; i < n - 1; ++i)
      if (i < irevs[i])
	std::swap(a[i], a[irevs[i]]);
  };

  void idft(Complex* a) override {
    auto iw = ws.begin();
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int j = 0; j < l; ++j) {
	Complex w = (*iw++).conj();
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

    for (int i = 1; i < n - 1; ++i)
      if (i < irevs[i])
	std::swap(a[i], a[irevs[i]]);
    double inv = 1.0 / n;
    for (int i = 0; i < n; ++i)
      a[i] *= inv;
  };

  std::vector<int> irevs;
  std::vector<Complex> ws;
};
