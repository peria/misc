#pragma once

#include <cmath>
#include <vector>

#include "fft.h"

class StockhamDIT final : public FFT {
 public:
  StockhamDIT() = default;
  ~StockhamDIT() override = default;

  const char* name() const override { return "StockhamDIT"; }

  void setUp(int n_) override {
    FFT::setUp(n_);

    work.resize(n);
    const double theta = -2 * M_PI / n;
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      const double t0 = theta * l;
      for (int k = 0; k < m; ++k) {
	double t = t0 * k;
	Complex w {std::cos(t), std::sin(t)};
	ws.push_back(w);
      }
    }
  }

  void tearDown() override {
    ws.clear();
    work.clear();
  }

  void dft(Complex* a) override {
    Complex* x = a;
    Complex* y = work.data();
    auto iw = ws.begin();
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int k = 0; k < m; ++k) {
	Complex w = *iw++;
	for (int j = 0; j < l; ++j) {
	  int ix0 = k * l + j;
	  int ix1 = ix0 + l * m;
	  int iy0 = k * 2 * l + j;
	  int iy1 = iy0 + l;
	  Complex x0 = x[ix0];
	  Complex x1 = x[ix1];
	  y[iy0] = x0 + x1;
	  y[iy1] = (x0 - x1) * w;
	}
      }
      std::swap(x, y);
    }
    if (x != a)
      for (int i = 0; i < n; ++i)
	a[i] = x[i];
  };

  void idft(Complex* a) override {
    Complex* x = a;
    Complex* y = work.data();
    auto iw = ws.begin();
    for (int l = 1, m = n / 2; l < n; l *= 2, m /= 2) {
      for (int k = 0; k < m; ++k) {
	Complex w = (*iw++).conj();
	for (int j = 0; j < l; ++j) {
	  int ix0 = k * l + j;
	  int ix1 = ix0 + l * m;
	  int iy0 = k * 2 * l + j;
	  int iy1 = iy0 + l;
	  Complex x0 = x[ix0];
	  Complex x1 = x[ix1];
	  y[iy0] = x0 + x1;
	  y[iy1] = (x0 - x1) * w;
	}
      }
      std::swap(x, y);
    }
    if (x != a)
      for (int i = 0; i < n; ++i)
	a[i] = x[i];
    double inv = 1.0 / n;
    for (int i = 0; i < n; ++i)
      a[i] *= inv;
  };

  std::vector<Complex> ws;
  std::vector<Complex> work;
};
