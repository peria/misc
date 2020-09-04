#pragma once

#include "fft.h"

class StockhamDIF final : public FFT {
 public:
  StockhamDIF() = default;
  ~StockhamDIF() override = default;
  const char* name() const override { return "StockhamDIF"; }

  void setUp(int n_) override {
    n = n_;
    work.resize(n);
    const double theta = -2 * M_PI / n;
    for (int l = n / 2, m = 1; m < n; l /= 2, m *= 2) {
      const double t0 = theta * l;
      for (int k = 0; k < m; ++k) {
	double t = t0 * k;
	Complex w {std::cos(t), std::sin(t)};
	ws.push_back(w);
      }
    }
  };

  void tearDown() override {
    work.clear();
    ws.clear();
  }

  void dft(Complex* a) override {
    Complex* x = a;
    Complex* y = work.data();

    auto iw = ws.begin();
    for (int l = n / 2, m = 1; m < n; l /= 2, m *= 2) {
      for (int k = 0; k < m; ++k) {
	Complex w = *iw++;
	for (int j = 0; j < l; ++j) {
	  int ix0 = k * 2 * l + j;
	  int ix1 = ix0 + l;
	  int iy0 = k * l + j;
	  int iy1 = iy0 + l * m;
	  Complex x0 = x[ix0];
	  Complex x1 = x[ix1] * w;
	  y[iy0] = x0 + x1;
	  y[iy1] = x0 - x1;
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
    for (int l = n / 2, m = 1; m < n; l /= 2, m *= 2) {
      for (int k = 0; k < m; ++k) {
	Complex w = (*iw++).conj();
	for (int j = 0; j < l; ++j) {
	  int ix0 = k * 2 * l + j;
	  int ix1 = ix0 + l;
	  int iy0 = k * l + j;
	  int iy1 = iy0 + l * m;
	  Complex x0 = x[ix0];
	  Complex x1 = x[ix1] * w;
	  y[iy0] = x0 + x1;
	  y[iy1] = x0 - x1;
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

  int n;
  std::vector<Complex> ws;
  std::vector<Complex> work;
};
