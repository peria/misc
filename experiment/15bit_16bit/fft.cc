#include "fft.h"

#include <cassert>
#include <cmath>
#include <iostream>

Fft::Fft(const int64 n,
         const int64 log2n,
         const Radix radix,
         std::vector<double>& work,
         std::vector<double>& table)
    : n(n), log2n(log2n), radix(radix), work_(work), table_(table) {
  int64 width = 1;
  auto itr = table_.begin();
  for (int64 k = 0; k < log2n; ++k) {
    for (int64 i = 0; i < width; ++i) {
      double wt = M_PI * i / width;
      *itr++ = std::cos(wt);
      *itr++ = -std::sin(wt);
    }
    width *= 2;
  }
}

void Fft::run(const Direction dir, std::vector<double>& data) const {
  if (dir == Direction::Backward) {
    for (int64 i = 0; i < n; ++i) {
      data[2 * i + 1] = -data[2 * i + 1];
    }
  }

  // <Core>
  double* y = work_.data();
  double* x = data.data();
  double* table = table_.data();
  int width = 1, height = n;
  for (int i = 0; i < log2n; ++i) {
    height /= 2;
    if (i % 2) {
      Radix2(width, height, table, y, x);
    } else {
      Radix2(width, height, table, x, (i == log2n - 1) ? x : y);
    }
    width *= 2;
    table += width;
  }
  if (radix == Radix::Three) {
    height /= 3;
    Radix3(width, height, table, x, x);
    width *= 3;
    table += width;
  }
  if (radix == Radix::Five) {
    height /= 5;
    Radix5(width, height, table, x, x);
    width *= 5;
    table += width;
  }
  // </Core>

  if (dir == Direction::Backward) {
    double inv = 1.0 / n;
    for (int64 i = 0; i < n; ++i) {
      data[2 * i] *= inv;
      data[2 * i + 1] *= -inv;
    }
  }
}

void Fft::Radix2(const int width,
                 const int height,
                 double* ptr,
                 double* x,
                 double* y) const {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;
      double wr = ptr[2 * i], wi = ptr[2 * i + 1];
      double tr = x[2 * ix1] * wr - x[2 * ix1 + 1] * wi;
      double ti = x[2 * ix1] * wi + x[2 * ix1 + 1] * wr;
      y[2 * iy1] = x[2 * ix0] - tr;
      y[2 * iy1 + 1] = x[2 * ix0 + 1] - ti;
      y[2 * iy0] = x[2 * ix0] + tr;
      y[2 * iy0 + 1] = x[2 * ix0 + 1] + ti;
    }
  }
}

void Fft::Radix3(const int width,
                 const int height,
                 double* ptr,
                 double* x,
                 double* y) const {
  // TODO: Implement
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;
      double wr = ptr[2 * i], wi = ptr[2 * i + 1];
      double tr = x[2 * ix1] * wr - x[2 * ix1 + 1] * wi;
      double ti = x[2 * ix1] * wi + x[2 * ix1 + 1] * wr;
      y[2 * iy1] = x[2 * ix0] - tr;
      y[2 * iy1 + 1] = x[2 * ix0 + 1] - ti;
      y[2 * iy0] = x[2 * ix0] + tr;
      y[2 * iy0 + 1] = x[2 * ix0 + 1] + ti;
    }
  }
}

void Fft::Radix5(const int width,
                 const int height,
                 double* ptr,
                 double* x,
                 double* y) const {
  // TODO: Implement
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;
      double wr = ptr[2 * i], wi = ptr[2 * i + 1];
      double tr = x[2 * ix1] * wr - x[2 * ix1 + 1] * wi;
      double ti = x[2 * ix1] * wi + x[2 * ix1 + 1] * wr;
      y[2 * iy1] = x[2 * ix0] - tr;
      y[2 * iy1 + 1] = x[2 * ix0 + 1] - ti;
      y[2 * iy0] = x[2 * ix0] + tr;
      y[2 * iy0 + 1] = x[2 * ix0 + 1] + ti;
    }
  }
}

void Rft::run(const Direction dir, std::vector<double>& x) const {
  if (dir == Direction::Backward) {
    // Convert from real to complex
    double t = 2 * M_PI / n;
    for (int k = 1, j = n / 2 - 1; k < n / 4; ++k, --j) {
      double x0r = x[2 * k], x0i = x[2 * k + 1];
      double x1r = x[2 * j], x1i = x[2 * j + 1];
      double ar = x0r - x1r, ai = x0i + x1i;
      double wr = 1 + sin(k * t);
      double wi = -cos(k * t);
      double br = (ar * wr - ai * wi) * 0.5;
      double bi = (ar * wi + ai * wr) * 0.5;
      x[2 * k] = x0r - br;
      x[2 * k + 1] = x0i - bi;
      x[2 * j] = x1r + br;
      x[2 * j + 1] = -(-x1i + bi);
    }
    double xr = x[0], xi = x[1];
    x[0] = (xr + xi) * 0.5;
    x[1] = (xr - xi) * 0.5;
  }

  Fft::run(dir, x);

  if (dir == Direction::Forward) {
    // Convert from complex to real
    double t = 2 * M_PI / n;
    for (int k = 1, j = n / 2 - 1; k < n / 4; ++k, --j) {
      double x0r = x[2 * k], x0i = x[2 * k + 1];
      double x1r = x[2 * j], x1i = x[2 * j + 1];
      double ar = x0r - x1r, ai = x0i + x1i;
      double wr = 1 + std::sin(k * t);
      double wi = std::cos(k * t);
      double br = (ar * wr - ai * wi) * 0.5;
      double bi = (ar * wi + ai * wr) * 0.5;
      x[2 * k] = x0r - br;
      x[2 * k + 1] = x0i - bi;
      x[2 * j] = x1r + br;
      x[2 * j + 1] = -(-x1i + bi);
    }
    double xr = x[0], xi = x[1];
    x[0] = xr + xi;
    x[1] = xr - xi;
  }
}

std::ostream& operator<<(std::ostream& ost, const Fft::Radix& r) {
  switch (r) {
  case Fft::Radix::Two:
    ost << "2";
    break;
  case Fft::Radix::Three:
    ost << "3";
    break;
  case Fft::Radix::Five:
    ost << "5";
    break;
  }
  return ost;
}
