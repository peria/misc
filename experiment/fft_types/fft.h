#pragma once

#include <cmath>
#include "complex.h"

class FFT {
 public:
  virtual ~FFT() {}
  virtual const char* name() const = 0;
  virtual void setUp(int n_) { n = n_; }
  virtual void tearDown() {}
  virtual void dft(Complex*) = 0;
  virtual void idft(Complex*) = 0;
  virtual void rft(double* x) {
    dft(reinterpret_cast<Complex*>(x));
    double x0 = x[0];
    double x1 = x[1];
    x[0] = x0 + x1;
    x[1] = x0 - x1;
    double th = M_PI / n;
    for (int i = 1, j = n - 1; i < j; ++i, --j) {
      double xir = x[2 * i], xii = x[2 * i + 1];
      double xjr = x[2 * j], xji = x[2 * j + 1];
      double wr = 1 + std::sin(th * i), wi = std::cos(th * i);
      double yr = xir - xjr, yi = xii + xji;
      double zr = (wr * yr - wi * yi) * 0.5;
      double zi = (wr * yi + wi * yr) * 0.5;
      x[2 * i] = xir - zr;
      x[2 * i + 1] = xii - zi;
      x[2 * j] = xjr + zr;
      x[2 * j + 1] = -(-xji + zi);
    }
    x[n + 1] = -x[n + 1];
  }
  virtual void irft(double* x) {
    double x0 = x[0], x1 = x[1];
    x[0] = (x0 + x1) * 0.5;
    x[1] = (x0 - x1) * 0.5;
    double th = M_PI / n;
    for (int i = 1, j = n - 1; i < j; ++i, --j) {
      double xir = x[2 * i], xii = x[2 * i + 1];
      double xjr = x[2 * j], xji = x[2 * j + 1];
      double wr = 1 + std::sin(th * i), wi = -std::cos(th * i);
      double yr = xir - xjr, yi = xii + xji;
      double zr = (wr * yr - wi * yi) * 0.5;
      double zi = (wr * yi + wi * yr) * 0.5;
      x[2 * i] = xir - zr;
      x[2 * i + 1] = xii - zi;
      x[2 * j] = xjr + zr;
      x[2 * j + 1] = xji - zi;
    }
    x[n + 1] = -x[n + 1];
    idft(reinterpret_cast<Complex*>(x));
  }

  int n;
};
