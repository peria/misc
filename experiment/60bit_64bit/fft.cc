#include "fft.h"

#include <cassert>
#include <cmath>
#include <iostream>

Fft::Fft(const int64 n, const int64 log2n, const Radix radix,
         std::vector<double> &work, std::vector<double> &table)
    : n(n), log2n(log2n), log4n(log2n > 1 ? 2 - (log2n + 2) % 3 : 0),
      log8n(log2n > 1 ? (log2n - log4n) / 3 : 0), radix(radix), work_(work),
      table_(table) {
  int64 height = n;
  Complex *tbl = reinterpret_cast<Complex *>(table_.data());
  Complex *tbl0 = tbl;
  for (int64 i = 0; i < log8n; ++i) {
    height /= 8;
    setTable(8, height, tbl);
    tbl += 7 * height;
  }
  for (int64 i = 0; i < log4n; ++i) {
    height /= 4;
    setTable(4, height, tbl);
    tbl += 3 * height;
  }
  if (log2n == 1) {
    height /= 2;
    setTable(2, height, tbl);
  }
}

void Fft::run(const Direction dir, std::vector<double> &data) const {
  if (dir == Direction::Backward) {
    for (int64 i = 0; i < n; ++i) {
      data[2 * i + 1] = -data[2 * i + 1];
    }
  }

  // <Core>
  Complex *y = reinterpret_cast<Complex *>(work_.data());
  Complex *x = reinterpret_cast<Complex *>(data.data());
  Complex *table = reinterpret_cast<Complex *>(table_.data());
  bool data_in_x = true;
  int width = 1, height = n;
  for (int64 i = 0; i < log8n; ++i) {
    height /= 8;
    if (height > 1) {
      if (i % 2) {
        radix8(width, height, table, y, x);
      } else {
        radix8(width, height, table, x, y);
      }
      data_in_x = !data_in_x;
    } else {
      radix8(width, height, table, data_in_x ? x : y, x);
    }
    width *= 8;
    table += 7 * height;
  }
  for (int64 i = 0; i < log4n; ++i) {
    height /= 4;
    if (height > 1) {
      if (i % 2) {
        radix4(width, height, table, y, x);
      } else {
        radix4(width, height, table, x, y);
      }
      data_in_x = !data_in_x;
    } else {
      radix4(width, height, table, data_in_x ? x : y, x);
    }
    width *= 4;
    table += 3 * height;
  }
  if (log2n == 1) {
    height /= 2;
    if (n == 2) {
      radix2(height, table, x, x);
    } else {
      radix2(height, table, x, y);
      data_in_x = !data_in_x;
    }
    width *= 2;
  }

  if (radix == Radix::Five) {
    height /= 5;
    radix5(width, height, data_in_x ? x : y, x);
  }
  if (radix == Radix::Three) {
    height /= 3;
    radix3(width, height, data_in_x ? x : y, x);
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

void Fft::setTable(int64 r, const int64 height, Complex *table) {
  const double theta = -2.0 * M_PI / (r * height);
  for (int64 i = 0; i < height; ++i) {
    for (int64 j = 1; j < r; ++j) {
      double t = theta * i * j;
      table[i * (r - 1) + j - 1] = Complex{std::cos(t), std::sin(t)};
    }
  }
}

void Fft::radix2(const int height, Complex *ptr, Complex *x, Complex *y) const {
#define X(A, B) x[(A)*height + (B)]
#define Y(A, B) y[(A)*2 + (B)]
  Complex c0 = X(0, 0);
  Complex c1 = X(1, 0);
  Complex d0 = c0 + c1;
  Complex d1 = c0 - c1;
  Y(0, 1) = d1;
  Y(0, 0) = d0;
  for (int64 j = 1; j < height; ++j) {
    Complex w = ptr[j - 1];
    Complex c0 = X(0, j);
    Complex c1 = X(1, j);
    Complex d0 = c0 + c1;
    Complex d1 = c0 - c1;
    Y(j, 0) = d0;
    Y(j, 1) = w * d1;
  }
#undef X
#undef Y
}

void Fft::radix4(const int width, const int height, Complex *ptr, Complex *x,
                 Complex *y) const {
#define X(A, B, C) x[((A)*height + (B)) * width + (C)]
#define Y(A, B, C) y[((A)*4 + (B)) * width + (C)]
  for (int64 i = 0; i < width; ++i) {
    Complex c0 = X(0, 0, i);
    Complex c1 = X(1, 0, i);
    Complex c2 = X(2, 0, i);
    Complex c3 = X(3, 0, i);
    Complex d0 = c0 + c2;
    Complex d1 = c0 - c2;
    Complex d2 = c1 + c3;
    Complex d3 = (c1 - c3).i();
    Y(0, 3, i) = d1 - d3;
    Y(0, 2, i) = d0 - d2;
    Y(0, 1, i) = d1 + d3;
    Y(0, 0, i) = d0 + d2;
  }
  for (int64 j = 1; j < height; ++j) {
    Complex w1 = ptr[3 * j];
    Complex w2 = ptr[3 * j + 1];
    Complex w3 = ptr[3 * j + 2];
    for (int64 i = 0; i < width; ++i) {
      Complex c0 = X(0, j, i);
      Complex c1 = X(1, j, i);
      Complex c2 = X(2, j, i);
      Complex c3 = X(3, j, i);
      Complex d0 = c0 + c2;
      Complex d1 = c0 - c2;
      Complex d2 = c1 + c3;
      Complex d3 = (c1 - c3).i();
      Y(j, 0, i) = d0 + d2;
      Y(j, 1, i) = w1 * (d1 + d3);
      Y(j, 2, i) = w2 * (d0 - d2);
      Y(j, 3, i) = w3 * (d1 - d3);
    }
  }
#undef X
#undef Y
}

void Fft::radix8(const int width, const int height, Complex *ptr, Complex *x,
                 Complex *y) const {
#define X(A, B, C) x[((A)*height + (B)) * width + (C)]
#define Y(A, B, C) y[((A)*8 + (B)) * width + (C)]
  static constexpr double kC81 = 0.70710678118654752;

  for (int64 i = 0; i < width; ++i) {
    Complex c0 = X(0, 0, i);
    Complex c1 = X(1, 0, i);
    Complex c2 = X(2, 0, i);
    Complex c3 = X(3, 0, i);
    Complex c4 = X(4, 0, i);
    Complex c5 = X(5, 0, i);
    Complex c6 = X(6, 0, i);
    Complex c7 = X(7, 0, i);
    Complex d0 = c0 + c4;
    Complex d1 = c0 - c4;
    Complex d2 = c2 + c6;
    Complex d3 = (c2 - c6).i();
    Complex d4 = c1 + c5;
    Complex d5 = c1 - c5;
    Complex d6 = c3 + c7;
    Complex d7 = c3 - c7;
    Complex e0 = d0 + d2;
    Complex e1 = d0 - d2;
    Complex e2 = d4 + d6;
    Complex e3 = (d4 - d6).i();
    Complex e4 = kC81 * (d5 - d7);
    Complex e5 = kC81 * (d5 + d7).i();
    Complex e6 = d1 + e4;
    Complex e7 = d1 - e4;
    Complex e8 = d3 + e5;
    Complex e9 = d3 - e5;
    Y(0, 0, i) = e0 + e2;
    Y(0, 1, i) = e6 + e8;
    Y(0, 2, i) = e1 + e3;
    Y(0, 3, i) = e7 - e9;
    Y(0, 4, i) = e0 - e2;
    Y(0, 5, i) = e7 + e9;
    Y(0, 6, i) = e1 - e3;
    Y(0, 7, i) = e6 - e8;
  }
  for (int64 j = 1; j < height; ++j) {
    Complex w1 = ptr[7 * j];
    Complex w2 = ptr[7 * j + 1];
    Complex w3 = ptr[7 * j + 2];
    Complex w4 = ptr[7 * j + 3];
    Complex w5 = ptr[7 * j + 4];
    Complex w6 = ptr[7 * j + 5];
    Complex w7 = ptr[7 * j + 6];
    for (int64 i = 0; i < width; ++i) {
      Complex c0 = X(0, j, i);
      Complex c1 = X(1, j, i);
      Complex c2 = X(2, j, i);
      Complex c3 = X(3, j, i);
      Complex c4 = X(4, j, i);
      Complex c5 = X(5, j, i);
      Complex c6 = X(6, j, i);
      Complex c7 = X(7, j, i);
      Complex d0 = c0 + c4;
      Complex d1 = c0 - c4;
      Complex d2 = c2 + c6;
      Complex d3 = (c2 - c6).i();
      Complex d4 = c1 + c5;
      Complex d5 = c1 - c5;
      Complex d6 = c3 + c7;
      Complex d7 = c3 - c7;
      Complex e0 = d0 + d2;
      Complex e1 = d0 - d2;
      Complex e2 = d4 + d6;
      Complex e3 = (d4 - d6).i();
      Complex e4 = kC81 * (d5 - d7);
      Complex e5 = kC81 * (d5 + d7).i();
      Complex e6 = d1 + e4;
      Complex e7 = d1 - e4;
      Complex e8 = d3 + e5;
      Complex e9 = d3 - e5;
      Y(j, 0, i) = e0 + e2;
      Y(j, 1, i) = w1 * (e6 + e8);
      Y(j, 2, i) = w2 * (e1 + e3);
      Y(j, 3, i) = w3 * (e7 - e9);
      Y(j, 4, i) = w4 * (e0 - e2);
      Y(j, 5, i) = w5 * (e7 + e9);
      Y(j, 6, i) = w6 * (e1 - e3);
      Y(j, 7, i) = w7 * (e6 - e8);
    }
  }
#undef X
#undef Y
}

void Fft::radix3(const int width, const int height, Complex *x,
                 Complex *y) const {
  static constexpr double kC31 = 0.86602540378443865;
  static constexpr double kC32 = 0.5;
#define X(A, C) x[(A)*height * width + (C)]
#define Y(B, C) y[(B)*width + (C)]

  for (int64 i = 0; i < width; ++i) {
    Complex c0 = X(0, i);
    Complex c1 = X(1, i);
    Complex c2 = X(2, i);
    Complex d0 = c1 + c2;
    Complex d1 = c0 - kC32 * d0;
    Complex d2 = kC31 * (c1 - c2).i();
    Y(0, i) = c0 + d0;
    Y(1, i) = d1 + d2;
    Y(2, i) = d1 - d2;
  }
#undef X
#undef Y
}

void Fft::radix5(const int width, const int height, Complex *x,
                 Complex *y) const {
#define X(A, C) x[(A)*height * width + (C)]
#define Y(B, C) y[(B)*width + (C)]
  static constexpr double kC51 = 0.95105651629515357;
  static constexpr double kC52 = 0.61803398874989485;
  static constexpr double kC53 = 0.55901699437494742;
  static constexpr double kC54 = 0.25;

  for (int64 i = 0; i < width; ++i) {
    Complex c0 = X(0, i);
    Complex c1 = X(1, i);
    Complex c2 = X(2, i);
    Complex c3 = X(3, i);
    Complex c4 = X(4, i);
    Complex d0 = c1 + c4;
    Complex d1 = c2 + c3;
    Complex d2 = kC51 * (c1 - c4);
    Complex d3 = kC51 * (c2 - c3);
    Complex d4 = d0 + d1;
    Complex d5 = kC53 * (d0 - d1);
    Complex d6 = c0 - kC54 * d4;
    Complex d7 = d6 + d5;
    Complex d8 = d6 - d5;
    Complex d9 = (d2 + kC52 * d3).i();
    Complex d10 = (kC52 * d2 - d3).i();
    Y(0, i) = c0 + d4;
    Y(1, i) = d7 + d9;
    Y(2, i) = d8 + d10;
    Y(3, i) = d8 - d10;
    Y(4, i) = d7 - d9;
  }
#undef X
#undef Y
}

void Rft::run(const Direction dir, std::vector<double> &x) const {
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

std::ostream &operator<<(std::ostream &ost, const Fft::Radix &r) {
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
